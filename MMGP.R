training_MM = function(
    data,
    prior_mean = NULL,
    ini_hp_0 = NULL,
    ini_hp_i = NULL,
    kern_0 = "SE",
    kern_i = "SE",
    common_hp = TRUE,
    grid_inputs = NULL,
    pen_diag = 1e-10,
    n_iter_max = 25,
    cv_threshold = 1e-3)
{
  list_probe = data$ProbeID %>% unique()
  
  n = 1
  t1 = Sys.time()
  
  floop1 = function(i){
    
    cat('ID n°', n, ':', i, '\n \n')
    n <<- n+1
    
    sub_db = data %>%
      filter(ProbeID == i) %>% 
      select(- ProbeID)
    
    train_magma(data = sub_db,
                prior_mean = prior_mean,
                ini_hp_0 = ini_hp_0,
                ini_hp_i = ini_hp_i,
                kern_0 = kern_0,
                kern_i = kern_i,
                common_hp = common_hp,
                grid_inputs = grid_inputs,
                pen_diag = pen_diag,
                n_iter_max = n_iter_max,
                cv_threshold = cv_threshold) %>% 
      return()
  }
  train = sapply(list_probe, floop1, simplify = FALSE, USE.NAMES = TRUE)
  
  t2 = Sys.time()
  
  train[['Training_time']] = t2-t1
  
  return(train)
}

pred_MM = function(data,
                   trained_model = NULL,
                   grid_inputs = NULL,
                   hp = NULL,
                   kern = "SE",
                   hyperpost = NULL,
                   get_hyperpost = FALSE,
                   get_full_cov = FALSE, 
                   plot = TRUE)
{
  new_ID = data$ProbeID %>% unique()
  
  pred_magma(data %>% select(- ProbeID),
             trained_model = trained_model[[new_ID]],
             grid_inputs = grid_inputs,
             hp = hp,
             kern = kern,
             hyperpost = hyperpost,
             get_hyperpost = get_hyperpost,
             get_full_cov = get_full_cov,
             plot = plot) %>% 
    return()
}

pred_MM_loop = function(db, trained_model, sep_train_test = 5, test_inputs = 6)
{
  db = db %>% mutate(Input = round(Input, 5))
  floop_j = function(j)
  {
    pred_dummy = pred_MM(db %>%
                           filter(ProbeID %in% j) %>%
                           filter(ID %in% db$ID[[1]]),
                         trained_model = trained_model,
                         grid_inputs = union(unique(db$Input), test_inputs),
                         get_hyperpost = TRUE,
                         plot = FALSE)
    
    hyperpost = pred_dummy$hyperpost
    
    floop_i = function(i)
    {
      cat("Indiv n°", i, "Probe n°", j, '\n \n')
      ID_pred = i
      ProbeID_pred = j
      
      db_test = db %>% 
        filter(ProbeID %in% ProbeID_pred) %>% 
        filter(ID %in% ID_pred) %>% 
        filter(Input < sep_train_test)
      
      test_point = db %>% 
        filter(ProbeID %in% ProbeID_pred) %>%
        filter(ID %in% ID_pred) %>% 
        filter(Input > sep_train_test)
      
      ## If testing data is missing, predict at 6 years
      if(nrow(test_point) == 0){
        inputs = test_inputs
      } else{
        inputs = test_point$Input
      }
      
      ## If individual data is missing, return prediction of the hyperposterior
      if(nrow(db_test) == 0){
        
        pred = hyperpost$pred %>% 
          filter(Input %in% inputs) %>% 
          select(- Reference) %>% 
          mutate(ID = i, ProbeID=j) %>% 
          relocate(Input, .after = Var)
        
      } else{ ## Just predict normally otherwise
        
        pred = pred_MM(db_test,
                       trained_model = trained_model,
                       grid_inputs = inputs, 
                       kern = 'SE + LIN',
                       hyperpost = hyperpost,
                       plot = FALSE) %>% 
          mutate('ID' = i, 'ProbeID' = j) 
      }
      
      pred %>% return()
    }
    db$ID %>%
      unique %>% 
      lapply(floop_i) %>% 
      bind_rows() %>% 
      return()
  }
  db$ProbeID %>%
    unique %>% 
    lapply(floop_j) %>%
    bind_rows %>% 
    return()
}

pred_age = function(pred, db_coef, clock_fct = identity,  nb_sample = 10000){
  intercept = db_coef %>% filter(ProbeID == '(Intercept)') %>% pull(Coef)
  
  floop_j = function(j){
    floop = function(i){
      cat("Indiv n°", j, "Probe n°", i, '\n \n')
      
      coef_i = db_coef %>%
        filter(ProbeID == i) %>% 
        pull(Coef)
      
      pred_i = pred %>% 
        filter(ID == j, ProbeID == i)
      
      rnorm(nb_sample, pred_i$Mean, sqrt(pred_i$Var)) * coef_i
    }
    samples = pred$ProbeID %>%
      unique() %>% 
      sapply(floop) %>% 
      rowSums() %>% 
      `+`(intercept) %>% 
      clock_fct()
    
    tibble(ID = j, Sample = samples) %>% 
      return()
  }
  pred$ID %>% 
    unique() %>% 
    lapply(floop_j) %>% 
    bind_rows() %>% 
    return()
}

draw_pred = function(pred, nb_sample = 10){
  floop = function(i){
    pred %>% 
      mutate(Sample = i, Pred_draw = rnorm(nrow(pred), Mean, sqrt(Var))) %>% 
      return()
  }
  1:nb_sample %>% 
    lapply(floop) %>% 
    bind_rows() %>% 
    return()
}

MSE = function(obs, pred)
{
  input = obs %>% pull(Input)
  value = obs %>% pull(Output)
  
  mix_pred = pred %>%
    filter(Input %in% input) %>%
    pull(Mean)
  
  (value - mix_pred)^2 %>%
    mean() %>%
    return()
}

WCIC = function(obs, pred, level)
{
  t = obs %>% pull(Input)
  value = obs %>% pull(Output)
  
  mean = pred %>% filter(Input %in% t) %>% pull(Mean)
  sd = pred %>% filter(Input %in% t) %>% pull(Var) %>% sqrt
  
  CI_inf = mean - qnorm(1 - level/2) * sd
  CI_sup = mean + qnorm(1 - level/2) * sd
  
  100 * ((CI_inf < value) & (value < CI_sup)) %>%
    mean %>%
    return()
}

eval = function(db, mod, name = 'Multi-mean')
{ 
  db = db %>% mutate(Input = round(Input, 4))
  floop_i = function(i){
    ## Compute the hyper-posterior at the adequate inputs
    pred_dummy = pred_MM(db %>%
                           filter(ProbeID %in% i) %>%
                           filter(ID %in% db$ID[[1]]),
                         trained_model = mod,
                         grid_inputs = unique(db$Input),
                         get_hyperpost = TRUE,
                         plot = FALSE)
    
    hyperpost = pred_dummy$hyperpost
    
    floop_j = function(j, hyp){
      cat("Probe n°", i, '- ID n°',j, '\n \n')
      
      ID_pred = j
      ProbeID_pred = i
      
      db_test = db %>% 
        filter(ProbeID %in% ProbeID_pred) %>% 
        filter(ID %in% ID_pred) %>% 
        filter(Input < 5) %>% 
        dplyr::select(- ID)
      
      test_point = db %>% 
        filter(ProbeID %in% ProbeID_pred) %>%
        filter(ID %in% ID_pred) %>% 
        filter(Input > 5)
      
      if((nrow(test_point) > 0) & (nrow(db_test) > 0 )){
        pred = pred_MM(db_test,
                       trained_model = mod,
                       grid_inputs = unique(test_point$Input),
                       hyperpost = hyperpost,
                       plot = FALSE)
        
        mse = test_point %>% MSE(pred)
        wcic = test_point %>% WCIC(pred, 0.05)
        
      } else if(nrow(test_point) > 0){
        
        pred = hyperpost$pred %>% filter(Input %in% test_point$Input)
        
        mse = test_point %>% MSE(pred)
        wcic = test_point %>% WCIC(pred, 0.05)
      } else {
        mse = NA
        wcic = NA
      }
      eval = tibble('Method' = name,
                    'ProbeID' = j, 
                    'MSE' = mse,
                    'WCIC' = wcic)
      eval %>% # bind_rows(eval_clust) %>%
        return()
    } 
    db$ID %>%
      unique %>% 
      lapply(floop_j, hyp = hyperpost) %>%
      bind_rows %>% 
      mutate(ID = i, .before = 3) %>% 
      return()
  } 
  db$ProbeID %>%
    unique %>% 
    lapply(floop_i) %>%
    bind_rows %>% 
    return()
}

reformat = function(db)
{
  db %>% 
    pivot_longer(cols = - ProbeID) %>% 
    mutate(ID = str_match(name, "-(.*?)$")[,2],
           Input = str_match(name, "M(.*?)-")[,2] %>% as.numeric(),
           Output = value) %>% 
    select(-c(name, value)) %>% 
    return()
}

anti.trafo = function(x, adult.age=20){
  
  ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) 
}
