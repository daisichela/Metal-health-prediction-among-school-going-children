pacman::p_load(tidyverse, gtsummary, magrittr, labelled, dplyr, mice,
               glmnet, tidymodels, ranger,xgboost, purrr,tibble,
               yardstick, tabpfn)
##Data Cleaning
{
##Reading in the data
data <- readRDS("~/PHD work/KIPMAT followup data/KIPMAT follow up/Data/Kipmat_KEMRI/combined_data_full.rds")

  ossdx <- data %>% 
    dplyr::select(
      dplyr::contains("ossdx"),
      dplyr::all_of(c(
        "study_number",
        "csex",
        "cage",
        "ulemavu",
        "specify_ulemavu",
        "kufitika",
        "relig",
        "edulev",
        "occup"
      ))
    )
  
  write.csv(ossdx, file = "ossdx all.csv")

  
  add_ossdx<- read.csv("~/PHD work/KIPMAT followup data/Data/UzimaKipmat_DATA_NOHDRS_2026-10-27_1510.csv")

#Visualizing missing data
visdat::vis_miss(data[,1:30])
visdat::vis_miss(data[,31:60])
visdat::vis_miss(data[,61:90])
visdat::vis_miss(data[,91:150])
visdat::vis_miss(data[,150:250])
visdat::vis_miss(data[,251:350])
visdat::vis_miss(data[,351:450])
visdat::vis_miss(data[,451:500])

data <- data %>% mutate(across(where(is.character), to_factor)) ## Converting character variables to factors
data<- val_labels_to_na(data)

##Function to rename  factor levels to NA when they are named as "Missing" or 9 
replace_with_na <- function(x, level_to_replace) {
  if (is.factor(x)) {
    levels(x)[levels(x) == level_to_replace] <- NA
  }
  return(x)
}
data <- as.data.frame(lapply(data, replace_with_na, level_to_replace = "Missing")) # Apply the function to each column in the data frame
data <- as.data.frame(lapply(data, replace_with_na, level_to_replace = "")) 

###calculaing maternal age

data$dobn <- as.Date(data$dobn, format = "%d-%b-%y")
data$dob <- as.Date(data$dob, format = "%d-%b-%y")
data %<>% mutate(mother_age = interval(data$dob, data$dobn)%/% months(1))#computing mother age from date of birth and date of baby's admission 
data$mother_age <- data$mother_age/12
data$mother_age[data$mother_age<=0]<-NA ##code as NA if age is 0


##Turning outliers into NA
plot(data$baby_weight)

data$baby_weight[data$baby_weight<100]<-NA
data$gestation [data$gestation>45]<-NA
data$gestation [data$gestation<20]<-NA
summary(data$baby_weight)

data %<>% mutate(occmothern = case_when(occmothern== "Farmer" ~ "Unemployed",
                                        occmothern== "Housewife" ~ "Unemployed",
                                        occmothern== "Student" ~ "Unemployed",
                                        occmothern== "Health sector" ~ "Employed",
                                        occmothern== "Other private sector employment" ~ "Employed",
                                        occmothern== "Public sector employment" ~ "Employed",
                                        occmothern== "Research sector" ~ "Employed",
                                        occmothern== "Business" ~ "Employed"))


data %<>% mutate(ancattend= case_when(ancattend==0~ "None",
                                                               ancattend==1 ~ "<4",
                                                               ancattend==2 ~ "<4",
                                                               ancattend==3 ~ "<4",
                                                               ancattend==4 ~ ">4",
                                                               ancattend==5 ~ ">4",
                                                               ancattend == ">5" ~ ">4"
))


data %<>% mutate(hiv_pregtest= case_when(hiv_pregtest== "ND" ~ NA,
                                         hiv_pregtest== "NR" ~ "Negative",
                                         hiv_pregtest== "R" ~ "Positive" ))

data %<>% mutate(delivemode= case_when(delivemode==1 ~ "Vaginal",
                                       delivemode==2 ~ "Vaginal",
                                       delivemode==4 ~ "CS",
                                       delivemode==5 ~ "CS",
                                       delivemode==12 ~ NA,
                                       delivemode== 14 ~ NA,
                                       delivemode== 15 ~ NA,
                                       delivemode==45 ~ NA
))

data %<>% mutate(childclass	= case_when(childclass	==0 ~ "Kindergarten",
                                        childclass	==1 ~ "Grade 1",
                                        childclass	==2 ~ "Grade 2",
                                        childclass	==3 ~ "Grade 3",
                                        childclass	==4 ~ "Grade 4",
                                        childclass	==5 ~ "Grade 5",
                                        childclass	== 6 ~ "Grade 6",
                                        childclass	== 7 ~ "Grade 7",
                                        childclass	==8 ~ "Grade 8",
))

data %<>% mutate(marst	= case_when(marst	==1 ~ "Not married",
                                   marst	==2 ~ "Married",
                                   marst	==3 ~ "Not married",
                                   marst	==4 ~ "Not married",
                                   marst	==5 ~ "Not married"
))


data %<>% mutate(occup	= case_when(occup	==1 ~ "Unemployed",
                                   occup	==2 ~ "Employed",
                                   occup	==3 ~ "Employed",
                                   occup	==4 ~ "Unemployed",
                                   occup	==5 ~ "Employed",
                                   occup ==6 ~ "Unemployed"
))


data %<>% mutate(edulev	= case_when(edulev	==1 ~ "None",
                                    edulev	==2 ~ "Primary",
                                    edulev	==3 ~ "Primary",
                                    edulev	==4 ~ "Secondary",
                                    edulev	==5 ~ "Secondary",
                                    edulev ==6 ~ "Tertiary",
                                    edulev ==7 ~ "Tertiary"
))
 
num<- c("cage","aheight", "aheight_2", "aheight_3", "aweight", "aweight_2",
        "aweight_3", "amuac", "amuac_2", "ahc", "ahc_2")


data [num] <- lapply(data[num], function(x) as.numeric(as.character(x)))

data [num] <- lapply(data[num], function(x) round(x,0))
#Extracting the year variable

data$year <- format(data$dobn, "%Y")
data$year <- as.numeric(format(data$dobn, "%Y"))
#selecting variables

data_child<- select(data,
                                ##birth outcomes
                               year,
                                baby_weight,
                                gestation,
                    
                                ##Demographics characteristics
                                occmothern,
                                educlevel,
                                ##Maternal history
                                mother_age,
                                firstpreg, 
                                npregs,
                                nsborn,
                                ancattend,
                                nmalaria,
                                ntt,
                                hbfirsttest,
                                hiv_pregtest,
                                supspec_iron,
                                supspec_folicacid,
                                ##pregnacy comlpications history
                                complhist_pet,
                                #Baby characteristics
                                sexn,
                                multiple, 
                                #Newborn and Delivery characteristics
                                newb_hr,
                               
                                apgar5,
                                
                                vitkgiven,
                                chlorhexidine,
                                nevirapine,
                                delivemode,
                                
                                cbcl_sch_1:total_cbcl,
                                comment_cbcl,
                                ravens_col_a1_response: child_has_autism ) 

data_child$nmalaria<- as.numeric((as.character(data_child$nmalaria)))
data_child$nsborn <- as.numeric(as.character(data_child$nsborn))


#Exploratory analysis
#Pregnancy to birth characteristics
data_child %>%
  select(c(baby_weight: delivemode)) %>% 
  tbl_summary(
    type = list(
      nsborn ~"continuous",
      nmalaria ~ "continuous",
      apgar5 ~ "continuous"
    ),
    label = list(
      baby_weight ~ "Birth weight",
      gestation ~ "Gestation at birth",
      occmothern ~ " Mother's occupation during pregnancy",
      educlevel ~ "Mother's education level during pregnancy",
      mother_age ~ "Age of the mother",
      firstpreg ~ "First born child", 
      firstpreg ~ "Number of previous pregnancies",
      nsborn ~ "Number of previous stillbirths",
      ancattend ~ "Antenatal care visits attended",
      nmalaria ~ " Number of malaria prophylaxis",
      ntt ~ "Tetanus vaccines received",
      hbfirsttest ~ "First haemoglobin test results during pregnancy",
      hiv_pregtest ~ "HIV status during pregnancy",
      hiv_pregtest ~ "Iron supplements during pregnancy",
      supspec_folicacid ~ "Folic acid supplements during pregnancy",
      complhist_pet ~ "History of high blood pressure",
      sexn ~ "Child sex",
      multiple ~ "Twin pregnancy",
      apgar5 ~ "Apgar score at 5 minutes",
      vitkgiven ~ "Vitamin K given at birth",
      chlorhexidine ~ "Chlorhexidine given at birth",
      nevirapine ~ "Nevirapine given",
      delivemode ~ "Delivery mode"
    )
  )




data_child$total_cbcl <- as.numeric(data_child$total_cbcl)
data_child %>%
  select(c(cbcl_sch_1: comment_cbcl)) %>% 
  tbl_summary()

#################CBCL coding###########
cbcl_missing <- data_child %>%
  select(contains("cbcl", ignore.case = TRUE)) %>%
  summarise(across(everything(), ~ mean(is.na(.)) * 100)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "missing_percent"
  ) %>%
  arrange(desc(missing_percent))

ggplot(cbcl_missing,
       aes(x = fct_reorder(variable, missing_percent),
           y = missing_percent)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "CBCL variable",
    y = "Missing data (%)",
    title = "Missingness across CBCL variables"
  ) +
  theme_minimal()  
## CBCL total problems by Group
##########################################################################
#CBCL subscale scores#
#########################

Anxious_items <- c("cbcl_sch_14","cbcl_sch_29","cbcl_sch_30","cbcl_sch_31",
                   "cbcl_sch_32","cbcl_sch_33","cbcl_sch_35","cbcl_sch_45",
                   "cbcl_sch_50","cbcl_sch_52","cbcl_sch_71","cbcl_sch_91","cbcl_sch_112")
Withdrawn_items <- c("cbcl_sch_5","cbcl_sch_42","cbcl_sch_65","cbcl_sch_69",
                     "cbcl_sch_75","cbcl_sch_102","cbcl_sch_103","cbcl_sch_111")
Somatic_items <- c("cbcl_sch_47","cbcl_sch_49","cbcl_sch_51","cbcl_sch_54",
                   "cbcl_sch_56a","cbcl_sch_56b","cbcl_sch_56c","cbcl_sch_56d",
                   "cbcl_sch_56e","cbcl_sch_56f","cbcl_sch_56g")
Social_items <- c("cbcl_sch_11","cbcl_sch_12","cbcl_sch_25","cbcl_sch_27",
                  "cbcl_sch_34","cbcl_sch_36","cbcl_sch_38","cbcl_sch_48",
                  "cbcl_sch_62","cbcl_sch_64","cbcl_sch_79")
Thought_items <- c("cbcl_sch_9","cbcl_sch_18","cbcl_sch_40","cbcl_sch_46","cbcl_sch_58","cbcl_sch_59","cbcl_sch_60",
                   "cbcl_sch_66","cbcl_sch_70","cbcl_sch_76","cbcl_sch_83","cbcl_sch_84","cbcl_sch_85","cbcl_sch_92","cbcl_sch_100")
Attention_items <- c(
  "cbcl_sch_1", "cbcl_sch_4", "cbcl_sch_8", "cbcl_sch_10", "cbcl_sch_13",
  "cbcl_sch_17", "cbcl_sch_41", "cbcl_sch_61", "cbcl_sch_78", "cbcl_sch_80"
)
Rule_items <- c("cbcl_sch_2","cbcl_sch_26","cbcl_sch_28","cbcl_sch_39","cbcl_sch_43",
                "cbcl_sch_63","cbcl_sch_67","cbcl_sch_72","cbcl_sch_73","cbcl_sch_81",
                "cbcl_sch_82","cbcl_sch_90","cbcl_sch_96","cbcl_sch_99",
                "cbcl_sch_101","cbcl_sch_105","cbcl_sch_106")
Aggressive_items <- c("cbcl_sch_3","cbcl_sch_16","cbcl_sch_19","cbcl_sch_20",
                      "cbcl_sch_21","cbcl_sch_22","cbcl_sch_23","cbcl_sch_37",
                      "cbcl_sch_57","cbcl_sch_68","cbcl_sch_86","cbcl_sch_87",
                      "cbcl_sch_88","cbcl_sch_89","cbcl_sch_94","cbcl_sch_95",
                      "cbcl_sch_97","cbcl_sch_104")
Other_items <- c("cbcl_sch_6","cbcl_sch_7","cbcl_sch_15","cbcl_sch_24","cbcl_sch_44",
                 "cbcl_sch_53","cbcl_sch_55","cbcl_sch_56h", "cbcl_sch_74","cbcl_sch_77","cbcl_sch_93",
                 "cbcl_sch_98","cbcl_sch_107","cbcl_sch_108","cbcl_sch_109","cbcl_sch_110")



## Calculate CBCL domain/subscale sum scores

data_child %<>%
  mutate(
    across(
      contains("cbcl_sch_"),
      ~ as.numeric(as.character(.))
    )
  )

data_child$cbcl_anxious <- rowSums(
  data_child[, Anxious_items],
  na.rm = FALSE
)

data_child$cbcl_withdrawn <- rowSums(
  data_child[, Withdrawn_items],
  na.rm = FALSE
)

data_child$cbcl_somatic <- rowSums(
  data_child[, Somatic_items],
  na.rm = FALSE
)

data_child$cbcl_social <- rowSums(
  data_child[, Social_items],
  na.rm = FALSE
)

data_child$cbcl_thought <- rowSums(
  data_child[, Thought_items],
  na.rm = FALSE
)

data_child$cbcl_attention <- rowSums(
  data_child[, Attention_items],
  na.rm = FALSE
)

data_child$cbcl_rule_breaking <- rowSums(
  data_child[, Rule_items],
  na.rm = FALSE
)

data_child$cbcl_aggressive <- rowSums(
  data_child[, Aggressive_items],
  na.rm = FALSE
)

data_child$cbcl_other <- rowSums(
  data_child[, Other_items],
  na.rm = FALSE
)
###Create total problems item list#
# Create Total Problems item list
Total_Problems_items <- unique(c(
  Anxious_items,
  Withdrawn_items,
  Somatic_items,
  Social_items,
  Thought_items,
  Attention_items,
  Rule_items,
  Aggressive_items,
  Other_items
))

length(Total_Problems_items)


##Compute total score#
# Total Problems score
data_child$cbcl_total_problems <- rowSums(
  data_child[, Total_Problems_items],
  na.rm = FALSE
)

summary(data_child$cbcl_total_problems)
 
cbcl_domains<- data_child %>% select(cbcl_anxious,cbcl_withdrawn,cbcl_somatic,
                      cbcl_social,cbcl_thought,cbcl_attention,cbcl_rule_breaking
                      ,cbcl_aggressive,cbcl_other)%>% 
  
  pivot_longer(
    cols = everything(),
    names_to = "Domain",
    values_to = "Score"
  )
ggplot(cbcl_domains, aes(x = Domain, y = Score, colour = Domain)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    width = 0.15,
    alpha = 0.30,
    size = 1
  ) +
  scale_x_discrete(labels= c("cbcl_anxious" = "Anxious","cbcl_withdrawn" = "Withdrawn/depressed","cbcl_somatic" ="Somatic complaints",
                             "cbcl_social" = "Social problems","cbcl_thought" ="Thought problems",
                             "cbcl_attention" = "Attention problems","cbcl_rule_breaking" ="Rule breaking behaviour"
                             ,"cbcl_aggressive" ="Agressive behaviour","cbcl_other" = "Other problems"))+
  scale_color_discrete(labels= c("cbcl_anxious" = "Anxious","cbcl_withdrawn" = "Withdrawn/depressed","cbcl_somatic" ="Somatic complaints",
                                 "cbcl_social" = "Social problems","cbcl_thought" ="Thought problems",
                                 "cbcl_attention" = "Attention problems","cbcl_rule_breaking" ="Rule breaking behaviour"
                                 ,"cbcl_aggressive" ="Agressive behaviour","cbcl_other" = "Other problems"))+
  labs(
    title = "",
    x = NULL,
    y = "CBCL score"
  ) +
  theme_minimal() +
  theme(
    strip.text           = element_text(face = "bold"),
    axis.text.x          = element_text(angle = 45, hjust = 1),
    legend.position      = "right",
    legend.justification = c("right", "top"),
    text                 = element_text(family = "Times New Roman", size = 15),
    panel.grid           = element_blank(),
    panel.grid.major     = element_blank(),
    panel.border         = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line            = element_line(color = "black"),
    axis.ticks           = element_line(color = "black")
  )

ggplot(data_child, aes(y = cbcl_total_problems, x="")) +
  geom_boxplot(outlier.shape = NA, color = "black",) +
  geom_jitter(
    width = 0.25,
    alpha = 0.50,
    size = 1,
    color = "maroon",
  ) +
  
  labs(
    title = "",
    x = NULL,
    y = "CBCL score"
  ) +
  theme_minimal() +
  theme(
    strip.text           = element_text(face = "bold"),
    text                 = element_text(family = "Times New Roman", size = 15),
    legend.position      = "none",
    panel.grid           = element_blank(),
    panel.grid.major     = element_blank(),
    panel.border         = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line            = element_line(color = "black"),
    axis.ticks           = element_line(color = "black")
  )

##Multiple imputation

data_child_imp<- select(data_child,
                    ##birth outcomes
                    year,
                    baby_weight,
                    gestation,
                    
                    ##Demographics characteristics
                    educlevel,
                    ##Maternal history
                    mother_age,
                    firstpreg, 
                    npregs,
                    ancattend,
                    nmalaria,
                    ntt,
                    hbfirsttest,
                    hiv_pregtest,
                    supspec_iron,
                    supspec_folicacid,
                    ##pregnacy comlpications history
                    complhist_pet,
                    #Baby characteristics
                    sexn,
                    multiple, 
                    #Newborn and Delivery characteristics
                    newb_hr,
                    
                    apgar5,
                    
                    vitkgiven,
                    nevirapine,
                    delivemode,
                    cbcl_anxious,cbcl_withdrawn,cbcl_somatic,
                    cbcl_social,cbcl_thought,cbcl_attention,cbcl_rule_breaking
                    ,cbcl_aggressive,cbcl_other, cbcl_total_problems) 

visdat::vis_miss(data_child_imp)
missing_data<- naniar::miss_var_summary(data_child_imp) 

#Removing entries which could not be traced in the original dataset

data_child_imp <- data_child_imp[-c(891:900),]
data_child_imp %<>%  mutate(across(where(is.character), as.factor)) 
data_child_imp$ancattend[data_child_imp$ancattend== "None"]<- NA
data_child_imp$ancattend <- droplevels(data_child_imp$ancattend)

#Splitting data before Multiple imputation

train_data <- data_child_imp[data_child_imp$year <=2018, ]
test_data  <- data_child_imp[data_child_imp$year == 2019, ]

#Number of observations in each dataset
n_train <- nrow(train_data)
n_test <- nrow(test_data)

#Setting an ignore argument that allows the imputation models 
#to be estimated using training 
#rows while still imputing missing values in the test rows
# FALSE = training rows used to estimate imputation models
# TRUE  = test rows excluded from model estimation
ignore_rows <- c(
  rep(FALSE, n_train),
  rep(TRUE, n_test)
)

# # --- Impute TRAINING set ---
init_all   <- mice(data_child_imp, maxit = 0)
pred_matrix  <- init_all$predictorMatrix

# # Exclude these variables from being used as predictors during imputation
pred_matrix[, c("cbcl_anxious","cbcl_withdrawn", "cbcl_somatic", "cbcl_social",
                "cbcl_thought", "cbcl_attention","cbcl_rule_breaking","cbcl_aggressive",
                "cbcl_other","cbcl_total_problems", "year")] <- 0

imp_all        <- mice(data_child_imp, maxit = 10, m = 5, predictorMatrix = pred_matrix,
                         ignore = ignore_rows,seed = 123)
imputed_all <- mice::complete(imp_all, action = "all")

imputed_train_list <- lapply(
  imputed_all,
  function(dat) dat[seq_len(n_train), , drop = FALSE]
)

imputed_test_list <- lapply(
  imputed_all,
  function(dat) dat[n_train + seq_len(n_test), , drop = FALSE]
)

# Save imputed training datasets
#purrr::walk2(
  #imputed_train_list,
 # seq_along(imputed_train_list),
 # ~ write.csv(.x, paste0("imputed_mental_health_training_data_", .y, ".csv"), row.names = FALSE)
#)

}

##Data import and model development
imputed_train_data<- lapply(1:5, function(i) read.csv(paste0("imputed_mental_health_training_data_", i,".csv"), stringsAsFactors = TRUE))

# Variables to remove
vars_to_remove <- c(
   "cbcl_anxious","cbcl_withdrawn",  "cbcl_social",
  "cbcl_thought", "cbcl_attention","cbcl_somatic", "cbcl_rule_breaking",
  "cbcl_aggressive", "cbcl_other"
)

#

 imputed_train_data <- lapply(imputed_train_data, function(df){
   df[,!names(df) %in% vars_to_remove, drop= FALSE]
 })

# Save imputed test datasets
#purrr::walk2(
 # imputed_test_list,
 # seq_along(imputed_test_list),
#  ~ write.csv(.x, paste0("imputed_mental_health_test_data_", .y, ".csv"), row.names = FALSE)
#)

imputed_test_data<- lapply(1:5, function(i) read.csv(paste0("imputed_mental_health_test_data_", i,".csv"), stringsAsFactors = TRUE))

imputed_test_data <- lapply(imputed_test_data, function(df){
  df[,!names(df) %in% vars_to_remove, drop= FALSE]
})
##Using LASOO regresiion for variable selection
{
train_imp1<- read.csv("imputed_mental_health_training_data_1.csv")
train_imp1$cbcl_total_problems <- scale(train_imp1$cbcl_total_problems)

outcome<- "cbcl_total_problems"
predictors <- c("baby_weight" ,
                "gestation",
                 "educlevel",
                "mother_age" ,
                "firstpreg" ,
                 "npregs",
                "ancattend",
                "nmalaria" ,
                 "ntt",
                 "hbfirsttest",
                 "hiv_pregtest",       
                "supspec_iron",
                 "supspec_folicacid",
                 "complhist_pet",      
                "sexn" ,
                 "multiple" ,
                "newb_hr" ,
                "apgar5",
                "vitkgiven" ,
                "nevirapine",
                "delivemode" )

lasso_fom <- reformulate(predictors)
x<- model.matrix(
  lasso_fom,
  data = train_imp1
) [,-1, drop = FALSE]


y<- train_imp1[[outcome]]

dim(x)
length(y)
# Fitting Lasso with ross validation

set.seed (123)

cv_lasso <- cv.glmnet(
  x=x,
  y=y,
  family ="gaussian",
  alpha=1,
  nfolds=5,
  standardize =
    TRUE,
  type.measure= "mse")
plot(cv_lasso)


#Extracting the selected variables

lasso_coeff <- coef(cv_lasso,
                    s="lambda.min")

lasso_coeff_df <- data.frame(
  term = rownames(lasso_coeff),
  coefficient= as.numeric(lasso_coeff)
)

selected_terms <- subset(lasso_coeff_df,
                         coefficient !=0 & term != "(Intercept)")




selected_original_vars <- predictors[
  vapply(predictors, function(v) {
    any(grepl(
      paste0("^", v),
      selected_terms$term
    ))
  }, logical(1))
]


variables_to_keep <- unique(c(
  outcome,
  selected_original_vars
))

imputed_train_data <- map(imputed_train_data,
                        ~select(.x, all_of(variables_to_keep)))

imputed_test_data <- map(imputed_test_data,
                          ~select(.x, all_of(variables_to_keep)))


names(imputed_train_data[[1]])
}

#Model development

#Linear regression
 lm_mod <- linear_reg() %>% 
   set_engine("lm")

#Elastic net regression

enet_mod <- linear_reg(
  penalty = tune(),
  mixture =tune()
) %>% 
  set_engine ("glmnet")

#Random forest

rf_mod <- rand_forest(
  mode = "regression",
  mtry = tune(),
  min_n = tune ()
) %>%
  
  set_engine(
    "ranger",
    importance = "impurity",
    splitrule = "variance",
    replace = FALSE
  )

#XGBoost
xgb_mod <- boost_tree(
  mode = "regression",
  trees = tune(),
  tree_depth = tune(),
  learn_rate = tune(),
  loss_reduction = tune(),
  sample_size = tune(),
  mtry = tune(),
  min_n = tune()
) %>% 
  set_engine ("xgboost")

#Function to fit models



run_regression_tuning_split <- function(
    train_data,
    test_data,
    model_name = c("lm","elastic_net", "rf", "xgb"),
    outcome = "cbcl_somatic",
    cv_folds = 5,
    n_bootstrap = 1000,
    ci_level = 0.95) {
  
  
  model_name <- match.arg (model_name)
  
  #Excluding children with missing outcome
  
  train_data <- train_data %>%
    dplyr::filter(!is.na(.data[[outcome]]))
  
  test_data <- test_data %>%
    dplyr::filter(!is.na(.data[[outcome]]))
  
  # ---------------------------------------------------------
  # Check the outcome and calculate scaling parameters
  # ---------------------------------------------------------
  
  if (!outcome %in% names(train_data) ||
      !outcome %in% names(test_data)) {
    stop("The outcome variable must exist in both datasets.")
  }
  
  if (anyNA(train_data[[outcome]]) ||
      anyNA(test_data[[outcome]])) {
    stop("The outcome contains missing values. Handle these before fitting.")
  }
  
  #outcome_mean <- mean(train_data[[outcome]])
  #outcome_sd   <- sd(train_data[[outcome]])
  
  #if (!is.finite(outcome_sd) || outcome_sd == 0) {
  #  stop("The training outcome must have a positive standard deviation.")
  #}
  

  # Scale using training data only
 # train_model <- train_data %>%
  #  mutate(
   #   cbcl_scaled = (.data[[outcome]] - outcome_mean) / outcome_sd
  #  ) %>%
  #  select(-all_of(outcome))
  
  # Keep the original test outcome for final evaluation
  #test_predictors <- test_data %>%
   # select(-all_of(outcome))


 
  train_model <- train_data 
  
  test_model <- test_data 
#Cross validation folds

set.seed(123)

train_folds <- vfold_cv(
  train_model,
  v= cv_folds
  )
  
#Preprocessing recipe

recipe_obj <- recipe(stats::reformulate(".", response = outcome),
                                        data = train_model) %>% 
  step_novel(all_nominal_predictors()) %>% 
  step_dummy(all_nominal_predictors())

#Select model and tuning paramaters

if (model_name == "lm") {
  
  mod <- parsnip::linear_reg() %>%
    parsnip::set_engine("lm")
  
} else if (model_name == "elastic_net") {
  mod <- enet_mod
  param_space <- parameters(
    penalty(),
    mixture()
  )
} else if (model_name == "rf"){
  mod <- rf_mod
  
  param_space <- parameters(
    mtry(),
    min_n()
  ) %>% 
    update (
      mtry = finalize (mtry (), train_model)
    )
  
} else {
  mod <- xgb_mod
  
  param_space <- extract_parameter_set_dials(mod) %>% 
    update(
      mtry = finalize (mtry(), train_model)
    )
}

  
##Biulding the workflow

workflow_obj <- workflow() %>% 
  add_model(mod) %>% 
  add_recipe (recipe_obj)

regression_metrics <- yardstick::metric_set(
  yardstick::rmse,
  yardstick:: mae,
  yardstick::rsq
)

if (model_name=="lm") {
  cv_res <- tune::fit_resamples(
    workflow_obj,
    resamples = train_folds,
    metrics = regression_metrics,
    control = tune::control_resamples(
      save_pred = TRUE
    )
  )
  
  best_params <- NULL
  
  final_fit <- workflows::fit(
    workflow_obj,
    data = train_model
  )
  
  cv_summary <- tune::collect_metrics(cv_res)
  
} else {
  
  set.seed(123)
  
  tune_res <- tune::tune_bayes(
    workflow_obj,
    resamples = train_folds,
    param_info = param_space,
    metrics = regression_metrics,
    initial = 5,
    iter = 20,
    control = tune::control_bayes(
      save_pred = TRUE,
      verbose = FALSE
    )
  )
  
  best_params <- tune::select_best(
    tune_res,
    metric = "rmse"
  )
  
  final_wf <- tune::finalize_workflow(
    workflow_obj,
    best_params
  )
  
  final_fit <- workflows::fit(
    final_wf,
    data = train_model
  )
  
  cv_summary <- tune::collect_metrics(tune_res)
}


#Bayesian hyperparameter tuning

set.seed(123)

tune_res <- tune_bayes (
  workflow_obj,
  resamples= train_folds,
  param_info = param_space,
  metrics = metric_set(rmse, mae, rsq_trad),
  initial =5,
  iter =20,
  control = control_bayes(
    save_pred = TRUE,
    verbose = FALSE
  )
)
  
#Selecting best parameters and final fit model

best_params <- select_best(
  tune_res,
  metric = "rmse"
  
)

final_wf <- finalize_workflow(
  workflow_obj,
  best_params
)

final_fit <- fit (
  final_wf,
  data= train_model
)

#Predict on the test dataset

test_predictors <- test_model %>%
  dplyr::select(-dplyr::all_of(outcome))

pred_original <- predict(
  final_fit,
  new_data = test_predictors
)$.pred

test_predictions <- tibble::tibble(
  truth = test_model[[outcome]],
  .pred = as.numeric(pred_original)
)

test_metrics <- regression_metrics(
  test_predictions,
  truth = truth,
  estimate = .pred
)

#Test set regression metrics

regression_metrics <- metric_set(
  rmse,
  mae,
  rsq_trad
)

test_metrics <- regression_metrics(
  test_predictions,
  truth = truth,
  estimate =.pred
)

#Bootstrap confidence intervals
set.seed(123)

boot_metrics <- map_dfr(seq_len(n_bootstrap), function(i){
  
 boot_sample <- test_predictions %>% 
   slice_sample (
     n=nrow (test_predictions),
     replace = TRUE
   )
 
 regression_metrics(
   boot_sample,
   truth = truth,
   estimate = .pred
 ) %>%
   select(.metric, .estimate) %>%
   mutate(iter = i)
})

alpha <- 1 - ci_level

ci_bounds <- boot_metrics %>%
  group_by(.metric) %>%
  summarise(
    lower = quantile(
      .estimate, alpha / 2, na.rm = TRUE
    ),
    upper = quantile(
      .estimate, 1 - alpha / 2, na.rm = TRUE
    ),
    .groups = "drop"
  )

test_metrics_with_ci <- test_metrics %>%
  left_join(ci_bounds, by = ".metric") %>%
  select(.metric, .estimate, lower, upper)

# ---------------------------------------------------------
# 11. Variable importance
# ---------------------------------------------------------

if (model_name == "rf") {
  
  rf_fit_obj <- extract_fit_parsnip(final_fit)$fit
  
  var_imp_df <- tibble(
    variable = names(rf_fit_obj$variable.importance),
    importance = unname(rf_fit_obj$variable.importance)
  ) %>%
    arrange(desc(importance))
  
} else if (model_name == "xgb") {
  
  xgb_fit_obj <- extract_fit_parsnip(final_fit)$fit
  
  var_imp_df <- xgboost::xgb.importance(
    model = xgb_fit_obj
  ) %>%
    as_tibble() %>%
    rename(
      variable = Feature,
      importance = Gain
    ) %>%
    arrange(desc(importance))
  
} else {
  
  # Elastic-net coefficients at the selected penalty
  glmnet_fit <- extract_fit_parsnip(final_fit)$fit
  
  coef_matrix <- as.matrix(
    coef(glmnet_fit, s = best_params$penalty)
  )
  
  var_imp_df <- tibble(
    variable = rownames(coef_matrix),
    coefficient = as.numeric(coef_matrix[, 1])
  ) %>%
    filter(variable != "(Intercept)", coefficient != 0) %>%
    mutate(importance = abs(coefficient)) %>%
    arrange(desc(importance))
}

# ---------------------------------------------------------
# 12. Return results
# ---------------------------------------------------------

list(
  model_name = model_name,
  final_fit = final_fit,
  best_params = best_params,
  cv_summary = collect_metrics(tune_res),
  test_metrics = test_metrics_with_ci,
  test_predictions = test_predictions,
  bootstrap_metrics = boot_metrics,
  variable_importance = var_imp_df
)
}
 



##Function to pool results


pool_prediction_results <- function(results, ci_level = 0.95) {
  
  library(dplyr)
  library(purrr)
  
  m <- length(results)
  
  if (m < 2) {
    stop("At least two imputation results are required.")
  }
  
  # Extract test metrics from each imputation
  performance_by_imp <- map2_dfr(
    results,
    seq_along(results),
    function(res, i) {
      res$test_metrics %>%
        select(.metric, .estimate) %>%
        mutate(imputation = i)
    }
  )
  
  # Descriptive summary across imputations
  performance_summary <- performance_by_imp %>%
    group_by(.metric) %>%
    summarise(
      n_imputations = sum(!is.na(.estimate)),
      mean_estimate = mean(.estimate, na.rm = TRUE),
      sd_between_imputations = sd(.estimate, na.rm = TRUE),
      min_estimate = min(.estimate, na.rm = TRUE),
      max_estimate = max(.estimate, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Extract bootstrap variance within each imputation
  pooling_data <- map2_dfr(
    results,
    seq_along(results),
    function(res, i) {
      
      estimates <- res$test_metrics %>%
        select(.metric, .estimate)
      
      within_variance <- res$bootstrap_metrics %>%
        group_by(.metric) %>%
        summarise(
          U = var(.estimate, na.rm = TRUE),
          .groups = "drop"
        )
      
      estimates %>%
        left_join(within_variance, by = ".metric") %>%
        mutate(imputation = i)
    }
  )
  
  # Approximate Rubin pooling
  alpha <- 1 - ci_level
  
  rubin_summary <- pooling_data %>%
    group_by(.metric) %>%
    summarise(
      m = sum(complete.cases(.estimate, U)),
      Q_bar = mean(.estimate[complete.cases(.estimate, U)]),
      U_bar = mean(U[complete.cases(.estimate, U)]),
      B = var(.estimate[complete.cases(.estimate, U)]),
      .groups = "drop"
    ) %>%
    mutate(
      total_variance = U_bar + (1 + 1 / m) * B,
      pooled_se = sqrt(total_variance),
      df = if_else(
        B > 0 & U_bar > 0,
        (m - 1) *
          (1 + U_bar / ((1 + 1 / m) * B))^2,
        Inf
      ),
      critical_value = if_else(
        is.finite(df),
        qt(1 - alpha / 2, df = df),
        qnorm(1 - alpha / 2)
      ),
      lower_ci = Q_bar - critical_value * pooled_se,
      upper_ci = Q_bar + critical_value * pooled_se
    ) %>%
    select(
      .metric, m, Q_bar, pooled_se,
      lower_ci, upper_ci, B, U_bar
    )
  
  list(
    performance_by_imputation = performance_by_imp,
    performance_summary = performance_summary,
    approximate_rubin_summary = rubin_summary
  )
}


# Running the models

# Elastic-net regression
results_enet <- lapply(seq_along(imputed_train_data),function(i){
  run_regression_tuning_split(
    train_data = imputed_train_data[[i]],
    test_data = imputed_test_data[[i]],
    model_name = "elastic_net",
    outcome = "cbcl_somatic"
  ) })


# Random forest regression
results_rf <- lapply(seq_along(imputed_train_data),function(i){
  run_regression_tuning_split(
    train_data = imputed_train_data[[i]],
    test_data = imputed_test_data[[i]],
    model_name = "rf",
    outcome = "cbcl_somatic"
  ) })


# XGBoost regression
results_xgb <-lapply(seq_along(imputed_train_data),function(i){
  run_regression_tuning_split(
    train_data = imputed_train_data[[i]],
    test_data = imputed_test_data[[i]],
    model_name = "xgb",
    outcome = "cbcl_somatic"
  ) })
# Elastic net
pooled_enet <- pool_prediction_results(results_enet)

# Random forest
pooled_rf <- pool_prediction_results(results_rf)

# XGBoost
pooled_xgb <- pool_prediction_results(results_xgb)



pooled_enet$performance_summary
pooled_rf$performance_summary
pooled_xgb$performance_summary


#Variable importance

# Elastic net: coefficient-based importance
results_enet[[1]]$variable_importance

# Random forest: model-based importance
results_rf[[1]]$variable_importance

# XGBoost: gain-based importance
results_xgb[[1]]$variable_importance

##Complete case analysis

train_data <- train_data %>% select(-any_of(vars_to_remove))
test_data <- test_data %>% select(-any_of(vars_to_remove))

train_data_complete<- train_data %>% drop_na()
test_data_complete <- test_data %>%  drop_na()

#Model fitting
results_enet_complete <- run_regression_tuning_split(
  train_data = train_data_complete,
  test_data  = test_data_complete,
  model_name = "elastic_net",
  outcome    = "cbcl_total_problems"
)
results_enet_complete$test_metrics



###Tabpfn code

nzchar(Sys.getenv("TABPFN_TOKEN"))
outcome <- "cbcl_total_problems"

predictors <- c(
  "baby_weight",
  "gestation",
  "educlevel",
  "mother_age",
  "firstpreg",
  "npregs" ,
  "ancattend",
  "nmalaria",
  "ntt",
  "hbfirsttest",
  "hiv_pregtest",
  "supspec_iron",
  "supspec_folicacid",
  "complhist_pet",
  "sexn",
  "multiple",
  "newb_hr",
  "apgar5",
   "vitkgiven" ,
  "nevirapine" ,
  "delivemode",
  "year"
)

predictors2 <- c("firstpreg",
                 "ancattend" ,
                 "nmalaria",
                 "hbfirsttest",
                 "complhist_pet",
                 "sexn" ,
                 "multiple",
                 "vitkgiven",
                 "year")
run_tabpfn_reg <- function(
    train_data,
    test_data,
    outcome = "cbcl_total_problems",
    predictors,
    device = "auto",
    seed = 123
){
  
  #Preparing predictors and outcomes
  
  x_train <- train_data[, predictors, drop=FALSE]
  x_test <- test_data[, predictors, drop= FALSE]
  
  y_train <- train_data[[outcome]]
  y_test <- test_data[[outcome]]
  
  #Fit model
  
  fitted_model <- tabpfn::tab_pfn(
    x_train,
    y_train,
    control = tabpfn::control_tab_pfn(
      device=device,
      random_state = seed,
      memory_saving_mode = "auto"
    )
  )
  
  #predict the held-out test data
  
  predictions <- predict(
    fitted_model,
    new_data= x_test)
  
  test_predictions <- tibble(
    truth = y_test,
    .pred = as.numeric(predictions$.pred)
  )
   #Evaluate on the original cbcl scale
  
  test_metrics <- yardstick::metric_set(
    yardstick::rmse,
    yardstick:: mae,
    yardstick::rsq
  )(
    test_predictions,
    truth = truth,
    estimate = .pred
  )
  list (
    model_name ="TabPFN",
    final_fit = fitted_model,
    test_predictions = test_predictions,
    test_metrics = test_metrics,
    n_train = nrow(x_train),
    n_test = nrow(x_test)
  )
}

##Run TabPFN model
results_tabpfn <- lapply(
  seq_along(imputed_train_data),
  function(i){
    message("Running TabPFN for imputation", i)
    
    run_tabpfn_reg(
      train_data= imputed_train_data[[i]],
      test_data = imputed_test_data[[i]],
      outcome = outcome,
      predictors = predictors2,
      device ="auto",
      seed = 123+i
    )
  }
    )

tabpfn_performance <- purrr::map2_dfr(
  results_tabpfn,
  seq_along(results_tabpfn),
  function(res, i) {
    res$test_metrics |>
      dplyr::mutate(imputation = i)
  }
)

tabpfn_summary <- tabpfn_performance |>
  dplyr::group_by(.metric) |>
  dplyr::summarise(
    mean_estimate = mean(.estimate, na.rm = TRUE),
    sd_between_imputations = sd(.estimate, na.rm = TRUE),
    min_estimate = min(.estimate, na.rm = TRUE),
    max_estimate = max(.estimate, na.rm = TRUE),
    .groups = "drop"
  )

tabpfn_summary


#Logistic regression

model <- lm(cbcl_total_problems~. ,data=imputed_train_data[[1]])
summary(model)

# Extract coefficients, confidence intervals and p-values
results <- tidy(model, conf.int = TRUE) %>%
  mutate(
    p_adjusted_BH = p.adjust(p.value, method = "BH"),
    significance = case_when(
      p_adjusted_BH < 0.001 ~ "***",
      p_adjusted_BH < 0.01  ~ "**",
      p_adjusted_BH < 0.05  ~ "*",
      TRUE ~ ""
    )
  ) %>%
  select(
    term, estimate, std.error, statistic,
    conf.low, conf.high, p.value,
    p_adjusted_BH, significance
  )

# Create a formatted table
results %>%
  gt() %>%
  cols_label(
    term = "Predictor",
    estimate = "Coefficient",
    std.error = "Standard error",
    statistic = "t-statistic",
    conf.low = "95% CI (lower)",
    conf.high = "95% CI (upper)",
    p.value = "Original p-value",
    p_adjusted_BH = "BH-adjusted p-value",
    significance = "Significance"
  ) %>%
  fmt_number(
    columns = c(
      estimate, std.error, statistic,
      conf.low, conf.high, p.value, p_adjusted_BH
    ),
    decimals = 3
  ) %>%
  tab_header(
    title = "Linear regression predicting CBCL total-problems scores",
    subtitle = "P-values adjusted using the Benjamini–Hochberg procedure"
  )


# Extract coefficient results
results <- as.data.frame(summary(model)$coefficients)

# Adjust p-values for multiple testing
results$p_adjusted_BH <- p.adjust(
  results$`Pr(>|t|)`,
  method = "BH"
)

# View results
results
tab_model(model)

