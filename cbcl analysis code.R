library(pacman)
p_load(epicalc,foreign,MASS,dplyr,readr,irr,readxl,MBESS,psych,ggplot2,leaps,haven,
       semTools,lavaan,semPlot,dplyr,nnet,ROCR,ggplot2,plotly,epiR,janitor,tidyr,stringr,Hmisc,multcomp)
setwd("C:/Users/ezra.too/OneDrive - Aga Khan University/Ezra/PhD/Data Analysis/ahos Data/Datasets")

ahos<-read_dta("C:/Ezra/PhD/Data Analysis/ahos Data/Datasets/ahos baseline dataset_21062022.dta")
ahos1<-read_dta("C:/Ezra/PhD/Data Analysis/ahos Data/Datasets/ahos PHASE1 Dataset.dta")
ahos2<-read_dta("C:/Ezra/PhD/Data Analysis/ahos Data/Datasets/ahos_PHASE2_Dataset.dta")
ahos_dat2 <- read.csv("ahos_baseline_cleaned_med_histparent_peer_rela.csv", header = T)
ahos_cd4_viral <- read_excel("CD4 Viral Load ahos Results.xlsx")
ahos_health_risk_behav <- read_excel("ACASI_cleaned.xlsx")

##Renaming part_id in ahos1 to studyno##
ahos<-ahos%>%
  rename(studyno = part_id)

# Merge "group" into wave 1
ahos1 <- ahos1 %>%
  left_join(ahos %>% select(studyno, group), by = "studyno")



#*            Descriptive statistics
#*            Age
mean(ahos1$cage)
sd(ahos1$cage)
summary(ahos1$cage)

age12<-ahos%>%filter(cage==12)
age11<-ahos%>%filter(cage==11)

age2<-ahos1%>%filter(cage==12)
age1<-ahos1%>%filter(cage==11)

#<><><><><><><><>Age
summary(ahos1$cage)

ahos1 %>%
  summarise(min_age = min(cage, na.rm =T),
            max_age = max(cage, na.rm =T),
            mean_age = mean(cage, na.rm =T),
            sd_age = sd(cage, na.rm =T))

ahos1 %>% group_by(group) %>%
  summarise(count = n(),mean=mean(cage), sd=sd(cage)) %>%
  mutate(summary = paste0(round(mean,2)," ","?"," ",round(sd,2)))

#Summary including age ranges across groups#
ahos1 %>%
  group_by(group) %>%
  summarise(
    count = n(),
    mean = mean(cage, na.rm = TRUE),
    sd = sd(cage, na.rm = TRUE),
    min_age = min(cage, na.rm = TRUE),
    max_age = max(cage, na.rm = TRUE),
    summary = paste0(
      round(mean(cage, na.rm = TRUE), 2),
      " ± ",
      round(sd(cage, na.rm = TRUE), 2)
    ),
    .groups = "drop"
  )


# Compute the analysis of variance
cage.aov <- aov(cage ~ group, data = ahos1)
summary(cage.aov)

#<><><><><><><><> HIV group

table(ahos1$group)

#<><><><><><><><>Education level 
table(ahos1$childclass)


library(dplyr)

ahos1 <- ahos1 %>%
  mutate(child_edu = case_when(
    childclass %in% c("", "N/A", "No","NA","None","Not at",
                      "Claris Masika","School drop uot","Lucy Mbura","Stephen Otienè onyango",
                      "Mother","Not at school","Lucy","NOT at school ") ~ "Not in school",
    childclass %in% c("Special class", "Special unit","Class 1 special unit ") ~ "lower-primmary (pre-primary to class 5)",
    childclass %in% c("1","4","53.1","Class 1","Class 2","Class 3","Class 4",
                      "class 5","Class 5","Class3","Class5","Grade 1","Kg 2",
                      "KG 2","Kg 3","PP 2","1year","Claas 5","PP1",
                      "Class2", "PP1","Class4","Class 1","class4,below averag(e",
                      "Class","Grade 4","Grade 5","C lass 5","Mwambani Primary School") ~ "lower-primmary (pre-primary to class 5)",
    childclass %in% c("Class 6","Class 7","Class 8","class 8","Class6","Class7",
                      "6","Class 8", "Class seven","Class8","7","Class Seven",
                      "C lass 6","class 7","CLASS 8","C lass 7") ~ "upper-primmary (class 6 to 8)",
    childclass %in% c("Form one","Form 1","Form 2","Form 3","Form 4",
                      "Form One", "Form1","Form2","Day Star Secondary School",
                      "Form 1M","Form onÃ","Form3") ~ "Secondary school",
    childclass %in% c("1st year","2nd year","Collage") ~ "Tertiary",
    TRUE ~ NA_character_
  ))

table(ahos1$child_edu)

#cross tabulation
ahos1 %>%
  group_by(group,child_edu) %>%
  summarise(n=n()) %>%
  mutate(prop=paste0("(",round(n/sum(n)*100,1),")"))

chisq.test(ahos1$child_edu,ahos1$group)# Chi square
fisher.test(ahos1$child_edu,ahos1$group, simulate.p.value = TRUE)

#<><><><><><><><>Parent alive (orphanhood) 
table(ahos1$parents_alive)

ahos1$parents_alive1 <- ifelse(ahos1$parents_alive == "A","Both parents are alive",
                              ifelse(ahos1$parents_alive == "B","Only the mother alive",
                                     ifelse(ahos1$parents_alive == "C","Only the father alive",
                                            ifelse(ahos1$parents_alive == "D","Both parents died","Both parents are alive"))))

table(ahos1$parents_alive1)

#cross tabulation
ahos1 %>%
  group_by(group,parents_alive1) %>%
  summarise(n=n()) %>%
  mutate(prop=paste0("(",round(n/sum(n)*100,1),")"))

chisq.test(ahos1$parents_alive1,ahos1$group)# Chi square

#<><><><><><><><>Sex 
table(ahos1$csex)

#cross tabulation
ahos1 %>%
  group_by(group,csex) %>%
  summarise(n=n()) %>%
  mutate(prop=paste0("(",round(n/sum(n)*100,1),")"))

chisq.test(ahos1$csex,ahos1$group)# Chi square



#<><><><><><><><>Social Economic Status 

summary(ahos1$SES)

ahos1 %>%
  summarise(min_ses = min(SES, na.rm =T),
            max_ses = max(SES, na.rm =T),
            mean_ses = mean(SES, na.rm =T),
            sd_ses = sd(SES, na.rm =T))


ahos1 %>% group_by(group) %>%
  summarise(count = n(),mean=mean(SES), sd=sd(SES),
            min_ses = min(SES, na.rm =T),
            max_ses = max(SES, na.rm =T)) %>%
  mutate(summary = paste0(round(mean,2)," ","?"," ",round(sd,2)))

#<><><><><><><><> Compute the analysis of variance
SES.aov <- aov(SES ~ group, data = ahos1)
summary(SES.aov)

#BMI for age  
#BMI = Weight (kg) / Height (m)?  

summary(ahos1[,c("aweight","aheight")]) # weight (kg), height (cm)

#Adding weight and height values for five ids

for(i in 1:dim(ahos1)[1]){
  if(is.na(ahos1$aweight[i])== T && ahos1$studyno[i] ==3147 ){ ahos1$aweight[i] <- 51.3}
  if(is.na(ahos1$aweight[i])== T && ahos1$studyno[i] ==3009 ){ ahos1$aweight[i] <- 50.9}
  if(is.na(ahos1$aweight[i])== T && ahos1$studyno[i] ==2141 ){ ahos1$aweight[i] <- 42.7}
  if(is.na(ahos1$aweight[i])== T && ahos1$studyno[i] ==2084 ){ ahos1$aweight[i] <- 36.6}
  if(is.na(ahos1$aweight[i])== T && ahos1$studyno[i] ==2106 ){ ahos1$aweight[i] <- 49.4}
  if(is.na(ahos1$aweight[i])== T ){ ahos1$aweight[i] <- 40.67} # replace missing values with mean
  
  if(is.na(ahos1$aheight[i])== T && ahos1$studyno[i] ==3147 ){ ahos1$aheight[i] <- 151}
  if(is.na(ahos1$aheight[i])== T && ahos1$studyno[i] ==3009 ){ ahos1$aheight[i] <- 168.6}
  if(is.na(ahos1$aheight[i])== T && ahos1$studyno[i] ==2141 ){ ahos1$aheight[i] <- 155.8}
  if(is.na(ahos1$aheight[i])== T && ahos1$studyno[i] ==2084 ){ ahos1$aheight[i] <- 147.1}
  if(is.na(ahos1$aheight[i])== T && ahos1$studyno[i] ==2106 ){ ahos1$aheight[i] <- 124}
  if(is.na(ahos1$aheight[i])== T ){ ahos1$aheight[i] <- 150.7} # replace missing values with mean
}



ahos1$aheight_m <- ahos1$aheight/100
ahos1$aweight1 <- ifelse(ahos1$aweight == 155.7, 55.7,ahos1$aweight)
summary(ahos1$aweight1)

summary(ahos1[,c("aweight","aheight_m")])
#ahos1[is.na(ahos1$aweight),"part_id"]

ahos1$BMI <- ahos1$aweight1/(ahos1$aheight_m^2)

summary(ahos1$BMI)
sd(ahos1$BMI)

#Summary including BMI ranges across groups#
ahos1 %>%
  group_by(group) %>%
  summarise(
    count = n(),
    mean = mean(BMI, na.rm = TRUE),
    sd = sd(BMI, na.rm = TRUE),
    min_BMI = min(BMI, na.rm = TRUE),
    max_BMI = max(BMI, na.rm = TRUE),
    summary = paste0(
      round(mean(BMI, na.rm = TRUE), 2),
      " ± ",
      round(sd(BMI, na.rm = TRUE), 2)
    ),
    .groups = "drop"
  )


# Compute the analysis of variance
BMI.aov <- aov(BMI ~ group, data = ahos1)
summary(BMI.aov)



###Weight classification#
#Underweight - BMI less than 18.5
#Normal healthy weight - BMI between 18.5 and 24.9
#Overweight - BMI between 25.0 and 29.9
#Obese - BMI between 30.0 and 39.9
#Morbidly obese - BMI 40.0 and above

ahos1$BMI_cat <- ifelse(ahos1$BMI < 18.5,"Underweight", 
                       ifelse(ahos1$BMI >= 18.5 & ahos1$BMI <= 24.9,"Normal healthy weigh","Overweight or Obese"))

table(ahos1$BMI_cat)


#cross tabulation
ahos1 %>%
  group_by(group,BMI_cat) %>%
  summarise(n=n()) %>%
  mutate(prop=paste0("(",round(n/sum(n)*100,1),")"))

chisq.test(ahos1$BMI_cat,ahos1$group)# Chi square
fisher.test(ahos1$BMI_cat,ahos1$group, simulate.p.value = TRUE)

##Stunting#
library(anthroplus)

ls("package:anthroplus")

# Convert age from years to months
ahos1$age_months <- ahos1$cage * 12

# Recode sex
ahos1$sex_code <- ifelse(ahos1$csex == "M", 1, 2)

# Calculate WHO 2007 height-for-age z-scores
haz <- anthroplus_zscores(
  sex = ahos1$sex_code,
  age = ahos1$age_months,
  height = ahos1$aheight,
  weight = ahos1$aweight1
)
names(haz)

##Attach HAZ to your dataset#
ahos1$HAZ <- haz$zhfa


# Define stunting
ahos1$stunted <- ifelse(
  is.na(ahos1$HAZ),
  "Missing",
  ifelse(ahos1$HAZ < -2, "Stunted", "Not stunted")
)

table(ahos1$stunted, useNA = "ifany")


##Alternative#
library(dplyr)

ahos1$stunted <- case_when(
  is.na(ahos1$HAZ) ~ "Missing",
  ahos1$HAZ < -2 ~ "Stunted",
  ahos1$HAZ >= -2 ~ "Not stunted"
)

table(ahos1$stunted, useNA = "ifany")

# Summary
table(ahos1$stunted)
prop.table(table(ahos1$stunted)) * 100


##Stunting by group#
ahos1 %>%
  group_by(group, stunted) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(percent = round(n / sum(n) * 100, 1))

tab_stunting <- table(ahos1$group, ahos1$stunted)

chisq.test(tab_stunting)

##Parenting#
summary(ahos1[,c("parenting1","parenting2","parenting3","parenting4","parenting5","parenting6","parenting7","parenting8","parenting9","parenting10",
                "parenting11","parenting12","parenting13","parenting14","parenting15","parenting16","parenting17","parenting18")])




ahos1$parenting <- apply(ahos1[,c("parenting1","parenting2","parenting3","parenting4","parenting5","parenting6","parenting7","parenting8","parenting9","parenting10",
                                "parenting11","parenting12","parenting13","parenting14","parenting15","parenting16","parenting17","parenting18")],1,sum ,na.rm = T)
summary(ahos1$parenting)
mean(ahos1$parenting)

##Group difference#
ahos1 %>%
  summarise(min_parenting = min(parenting, na.rm =T),
            max_parenting = max(parenting, na.rm =T),
            mean_parenting = mean(parenting, na.rm =T),
            sd_parenting = sd(parenting, na.rm =T))


ahos1 %>% 
  group_by(group) %>%
  summarise(
    count = n(),
    mean = mean(parenting, na.rm = TRUE),
    sd = sd(parenting, na.rm = TRUE),
    min_parenting = min(parenting, na.rm = TRUE),
    max_parenting = max(parenting, na.rm = TRUE)
  ) %>%
  mutate(summary = paste0(round(mean, 2), " ± ", round(sd, 2)))

parenting.aov <- aov(parenting ~ group, data = ahos1)
summary(parenting.aov)

##School belonging#
ahos1$sch_exp <- apply(ahos1[,c("sch_experience1","sch_experience2","sch_experience3","sch_experience4")],1,sum ,na.rm = T)
summary(ahos1$sch_exp)

##Group difference#
ahos1 %>%
  summarise(min_school = min(sch_exp, na.rm =T),
            max_school = max(sch_exp, na.rm =T),
            mean_school = mean(sch_exp, na.rm =T),
            sd_school = sd(sch_exp, na.rm =T))


ahos1 %>% 
  group_by(group) %>%
  summarise(
    count = n(),
    mean = mean(sch_exp, na.rm = TRUE),
    sd = sd(sch_exp, na.rm = TRUE),
    min_school = min(sch_exp, na.rm = TRUE),
    max_school = max(sch_exp, na.rm = TRUE)
  ) %>%
  mutate(summary = paste0(round(mean, 2), " ± ", round(sd, 2)))

school.aov <- aov(sch_exp ~ group, data = ahos1)
summary(school.aov)


##Peer relationships##
ahos1$peer_exp <- apply(ahos1[,c("peer_experience1","peer_experience2","peer_experience3","peer_experience4")],1,sum ,na.rm = T)
summary(ahos1$peer_exp)

# Reverse score item 1
ahos1$peer_experience1_rev <- 4 - ahos1$peer_experience1

# Compute summed score
ahos1$peer_exp <- apply(
  ahos1[, c("peer_experience1_rev",
           "peer_experience2",
           "peer_experience3",
           "peer_experience4")],
  1,
  sum,
  na.rm = TRUE
)

# Summary
summary(ahos1$peer_exp)

##Group difference#
ahos1 %>%
  summarise(min_peer = min(peer_exp, na.rm =T),
            max_peer = max(peer_exp, na.rm =T),
            mean_peer = mean(peer_exp, na.rm =T),
            sd_peer = sd(peer_exp, na.rm =T))


ahos1 %>% 
  group_by(group) %>%
  summarise(
    count = n(),
    mean = mean(peer_exp, na.rm = TRUE),
    sd = sd(peer_exp, na.rm = TRUE),
    min_school = min(peer_exp, na.rm = TRUE),
    max_school = max(peer_exp, na.rm = TRUE)
  ) %>%
  mutate(summary = paste0(round(mean, 2), " ± ", round(sd, 2)))

peer.aov <- aov(peer_exp ~ group, data = ahos1)
summary(peer.aov)


#<><><><><><><><>History of any other medical conditions
ahos1$medical_hist <- "no"
table(ahos1$he)	
table(ahos1$ccp)	
table(ahos1$ail)	
table(ahos1$scd)	
table(ahos1$stk)

ahos1$medical_history <- ifelse(ahos1$he == "y" | ahos1$ccp == "y" |  ahos1$ail == "y" |
                                 ahos1$scd == "y" |ahos1$stk == "y","yes", ahos1$medical_hist)

table(ahos1$medical_history)
table(ahos1$medical_hist)

##Past medical history by group#
ahos1 %>%
  group_by(group, medical_history) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(percent = round(n / sum(n) * 100, 1))

tab_hist <- table(ahos1$group, ahos1$medical_history)

chisq.test(tab_hist)

##Depression#
ahos1$phq_tt <- apply(ahos1[,c("phq1","phq2","phq3","phq4","phq5","phq6","phq7","phq8","phq9")],1,sum ,na.rm = T)
summary(ahos1$phq_tt)


ahos1$depression <- ifelse(ahos1$phq_tt >= 10,"Yes","No")
table (ahos1$depression)

##Depression by group#
ahos1 %>%
  group_by(group, depression) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(percent = round(n / sum(n) * 100, 1))

tab_dep <- table(ahos1$group, ahos1$depression)

chisq.test(tab_dep)

#<><><><><><><><>CD4 Cell count 
ahos_cd4_viral<-ahos_cd4_viral%>%
  rename(studyno = part_id)

ahos1 <- merge(ahos1,ahos_cd4_viral, by = "studyno", all.x = T)

summary(ahos1$CD4)

ahos1$CD4_cat <- ifelse(ahos1$CD4 < 350, "<350 cells/mm3", ifelse(ahos1$CD4 >= 350 & ahos1$CD4 <= 500, 
                                                               "350-500 cells/mm3","> 500 cells/mm3"))

ahos1 %>%
  group_by(group,CD4_cat) %>%
  summarise(n=n()) %>%
  mutate(prop=paste0("(",round(n/sum(n)*100,1),")"))


ahos1$CD4_cat1 <- ahos1$CD4_cat

for(i in 1:dim(ahos1)[1]){
  if(is.na(ahos1$CD4_cat[i]) ==T){ahos1$CD4_cat1[i] <- "None"}
}

table(ahos1$CD4_cat1)


#<><><><><><><><>HIV viral load concentration 
table(ahos1$Viral_Load)

regexp <- "[[:digit:]]+"
ahos1$Viral_Load1 = as.numeric(str_extract(ahos1$Viral_Load,regexp))

table(ahos1$Viral_Load1)
summary(ahos1$Viral_Load1)

ahos1$Viral_Load1_cat <- ifelse(ahos1$Viral_Load1 <= 1000,"< 1000 copies/ml",">1000 copies/ml")

table(ahos1$Viral_Load1_cat)



ahos1 %>%
  group_by(group,Viral_Load1_cat) %>%
  summarise(n=n()) %>%
  mutate(prop=paste0("(",round(n/sum(n)*100,1),")"))

ahos1$Viral_Load1_cat1 <- ahos1$Viral_Load1_cat

for(i in 1:dim(ahos1)[1]){
  if(is.na(ahos1$Viral_Load1_cat[i]) ==T){ahos1$Viral_Load1_cat1[i] <- "0None"}
}

#####################################################################
##Create suppressed, non-suprressed and invalid groups#
#######################################################################
library(dplyr)
library(stringr)
table(ahos1$Viral_Load)

ahos1 <- ahos1 %>%
  mutate(
    Viral_Load_clean = str_trim(Viral_Load),
    
    Viral_Load_numeric = as.numeric(
      str_extract(Viral_Load_clean, "\\d+\\.?\\d*")
    ),
    
    Viral_Category = case_when(
      
      # Missing values
      is.na(Viral_Load_clean) ~ "Missing",
      
      # Invalid results
      str_to_lower(Viral_Load_clean) == "invalid" ~ "Invalid",
      
      # Not detected / not detectable
      str_to_lower(Viral_Load_clean) %in%
        c("not detectable", "not detected") ~ "Suppressed",
      
      # Values reported as <x (e.g., <40)
      str_detect(Viral_Load_clean, "^<") ~ "Suppressed",
      
      # Numeric values <1000
      !is.na(Viral_Load_numeric) & Viral_Load_numeric < 1000 ~ "Suppressed",
      
      # Numeric values >=1000
      !is.na(Viral_Load_numeric) & Viral_Load_numeric >= 1000 ~ "Non-suppressed",
      
      TRUE ~ "Missing"
    )
  )

table(ahos1$Viral_Category, useNA = "ifany")

###Create a disease severity variable with controls as 0#
#0=control or exposed
#1=infected and virally suppressed
#2=infected and non-suppressed
#NA=infected but viral load missing or invalid

ahos1 <- ahos1 %>%
  mutate(
    disease_severity = case_when(
      
      # Non-infected groups
      group %in% c("control", "exposed") ~ 0,
      
      # Infected participants
      group == "infected" & Viral_Category == "Suppressed" ~ 1,
      group == "infected" & Viral_Category == "Non-suppressed" ~ 2,
      
      # Missing or invalid viral load among infected participants
      TRUE ~ NA_real_
    )
  )

ahos1$disease_severity <- factor(
  ahos1$disease_severity,
  levels = c(0, 1, 2),
  labels = c("Control/Exposed", "Suppressed", "Non-suppressed")
)

table(ahos1$disease_severity)

table(ahos1$group, ahos1$disease_severity, useNA = "ifany")



####################
#backward digit span #
#Total correct scores#
table(ahos1$digitspan_nrtc)

hist(ahos1$digitspan_nrtc)

summary(ahos1$digitspan_nrtc)
sd(ahos1$digitspan_nrtc)

table(ahos1$digitspan_nrtc, ahos1$group)

mean(ahos1$digitspan_nrtc == 0,na.rm = TRUE)

ahos1 %>%
  summarise(min_correct = min(digitspan_nrtc, na.rm =T),
            max_correct = max(digitspan_nrtc, na.rm =T),
            mean_correct = mean(digitspan_nrtc, na.rm =T),
            sd_correct = sd(digitspan_nrtc, na.rm =T))


ahos1 %>% 
  group_by(group) %>%
  summarise(
    count = n(),
    mean = mean(digitspan_nrtc, na.rm = TRUE),
    sd = sd(digitspan_nrtc, na.rm = TRUE),
    min_correct = min(digitspan_nrtc, na.rm = TRUE),
    max_correct = max(digitspan_nrtc, na.rm = TRUE)
  ) %>%
  mutate(summary = paste0(round(mean, 2), " ± ", round(sd, 2)))

#Highest level reached#
table(ahos1$digitspan_nrhlr)

barplot(table(ahos1$digitspan_nrhlr),
        xlab="Highest Level Reached",
        ylab="Frequency",
        main="Backward Digit Span Performance")

##Categorize highest level reached#
ahos1$digitspan_cat <- factor(
  ifelse(ahos1$digitspan_nrhlr < 3, "Low",
         ifelse(ahos1$digitspan_nrhlr <= 4, "Moderate", "High")),
  levels = c("Low", "Moderate", "High"),
  ordered = TRUE
)

table(ahos1$digitspan_cat)
prop.table(table(ahos1$digitspan_cat))

##Scores per group#
table(ahos1$digitspan_cat, ahos1$group)
tab <- table(ahos1$digitspan_cat, ahos1$group)
prop.table(tab, margin = 2)

round(prop.table(tab, margin = 2) * 100, 1)

barplot(table(ahos1$digitspan_cat),
        col = "lightblue",
        main = "Working Memory Categories",
        ylab = "Frequency")

chisq.test(table(ahos1$digitspan_cat, ahos1$group))


############################################
#Health Risk Behaviors#
##############################################
ahos_health_risk_behav <- read_excel("ACASI_cleaned.xlsx")

#Subsetting phase1 data from ahos1 HRBs data#
hrbs_phase1<-ahos_health_risk_behav[560:1098,]

#Subsetting phase2 data from ahos1 HRBs data#
hrbs_phase2<-ahos_health_risk_behav[1099:1558,]

#Merging with ahos1 data#
hrbs_phase1<-hrbs_phase1%>%
  rename(studyno = part_id)

#Remove 4 from study ids#
hrbs_phase1$studyno <- hrbs_phase1$studyno - 40000

#Verify that ids match#
range(ahos1$studyno, na.rm = TRUE)
range(hrbs_phase1$studyno, na.rm = TRUE)

#Check how many ids overlap#
length(intersect(ahos1$studyno, hrbs_phase1$studyno))

#Now merge#
ahos1 <- merge(ahos1,hrbs_phase1, by = "studyno", all.x = T)

#Sexual abuse#
#Lifetime sexual abuse#
table(ahos1$`_12`)

##Sexual abuse in past 12 months#
table(ahos1$`_14`)

##Composite sexual abuse score#
##Dichotomize item 14#
ahos1$`_14`[ahos1$`_14` == 999] <- NA
ahos1$`_12`[ahos1$`_12` == 999] <- NA

ahos1$q14_bin <- ifelse(ahos1$`_14` %in% c(1,2), 0,
                       ifelse(ahos1$`_14` %in% c(3,4,5,6), 1, NA))

table(ahos1$q14_bin)

#Recode item 12#
ahos1$q12_bin <- ifelse(ahos1$`_12` == 1, 1,
                       ifelse(ahos1$`_12` == 2, 0, NA))

##Composite sexual abuse variable#
ahos1$sexual_abuse_any <- pmax(ahos1$q12_bin, ahos1$q14_bin, na.rm = TRUE)
ahos1$sexual_abuse_any[is.infinite(ahos1$sexual_abuse_any)] <- NA

table(ahos1$sexual_abuse_any, useNA = "always")
prop.table(table(ahos1$sexual_abuse_any, useNA = "always"))

#Sexual abuse by group#
table(ahos1$sexual_abuse_any, ahos1$group, useNA = "always")
tab <- table(ahos1$sexual_abuse_any, ahos1$group, useNA = "always")
prop.table(tab, margin = 2)
round(prop.table(tab, margin = 2) * 100, 1)

chisq.test(table(ahos1$sexual_abuse_any, ahos1$group))


##Early sexual debut#
table(ahos1$`_38`)

ahos1$`_38`[ahos1$`_38` == 999] <- NA
ahos1$`_38`[ahos1$`_38` == 888] <- NA

#Categorize age of debut#
#0=never had sex, 1=Early debut,2=Later debut)#
ahos1$sex_debut_cat <- ifelse(ahos1$`_38` == 1, 0,
                             ifelse(ahos1$`_38` %in% 2:5, 1,
                                    ifelse(ahos1$`_38` %in% 6:9, 2, NA)))

ahos1$sex_debut_cat <- factor(ahos1$sex_debut_cat,
                             levels = c(0,1,2),
                             labels = c("Never had sex", "Early debut (≤14)", "Later debut (≥15)"))

table(ahos1$sex_debut_cat, useNA = "always")
prop.table(table(ahos1$sex_debut_cat,useNA = "always"))

#Sexual debut by group#
table(ahos1$sex_debut_cat, ahos1$group, useNA = "always")
tab <- table(ahos1$sex_debut_cat, ahos1$group, useNA = "always")
prop.table(tab, margin = 2)
round(prop.table(tab, margin = 2) * 100, 1)

chisq.test(table(ahos1$sex_debut_cat, ahos1$group))
fisher.test(table(ahos1$sex_debut_cat, ahos1$group))

##Transactional sex#
#Perpetrator#
ahos1$`_42`[ahos1$`_42` == 888] <- NA
table(ahos1$`_42`,useNA = "always")

#Victim#
ahos1$`_43`[ahos1$`_43` == 888] <- NA
ahos1$`_43`[ahos1$`_43` == 999] <- NA
table(ahos1$`_43`,useNA = "always")

#Composite transactional sex variable#
ahos1$tx_perp_bin <- ifelse(ahos1$`_42` == 1, 1,
                           ifelse(ahos1$`_42` == 2, 0, NA))

ahos1$tx_victim_bin <- ifelse(ahos1$`_43` == 1, 1,
                             ifelse(ahos1$`_43` == 2, 0, NA))


ahos1$transactional_sex_any <- pmax(ahos1$tx_perp_bin,
                                   ahos1$tx_victim_bin,
                                   na.rm = TRUE)

ahos1$transactional_sex_any[is.infinite(ahos1$transactional_sex_any)] <- NA

table(ahos1$transactional_sex_any, useNA = "always")
prop.table(table(ahos1$transactional_sex_any, useNA = "always"))

#Check overlap#
table(ahos1$tx_perp_bin, ahos1$tx_victim_bin)

#Transactional sex by group#
table(ahos1$transactional_sex_any, ahos1$group, useNA = "always")
tab <- table(ahos1$transactional_sex_any, ahos1$group, useNA = "always")
prop.table(tab, margin = 2)
round(prop.table(tab, margin = 2) * 100, 1)

chisq.test(table(ahos1$transactional_sex_any, ahos1$group))
fisher.test(table(ahos1$transactional_sex_any, ahos1$group))


#Tobacco use#
#Lifetime#
ahos1$`_21`[ahos1$`_21` == 999] <- NA
table(ahos1$`_21`,useNA = "always")

#Categorize lifetime tobacco use#
ahos1$tobacco_lifetime <- ifelse(ahos1$`_21` == 1, 0,
                                ifelse(ahos1$`_21` %in% 2:8 , 1, NA))

table(ahos1$tobacco_lifetime)

#Recent#
table(ahos1$`_22`,useNA = "always")

#Categorize recent tobacco use#
ahos1$tobacco_recent <- ifelse(ahos1$`_22` == 1, 0,
                              ifelse(ahos1$`_22` %in% 2:7, 1, NA))

table(ahos1$tobacco_recent)

#Alcohol use#
#Lifetime#
ahos1$`_25`[ahos1$`_25` == 999] <- NA
table(ahos1$`_25`,useNA = "always")

##Categorize lifetime alcohol use#
ahos1$alcohol_lifetime <- ifelse(ahos1$`_25` == 1, 0,
                                ifelse(ahos1$`_25` %in% 2:8 , 1, NA))

table(ahos1$alcohol_lifetime)

#Recent#
table(ahos1$`_26`,useNA = "always")

#Categorize recent alcohol use#
ahos1$alcohol_recent <- ifelse(ahos1$`_26` == 1, 0,
                              ifelse(ahos1$`_26` %in% 2:7, 1, NA))

table(ahos1$alcohol_recent)

##Composite alcohol/tobacco use#
ahos1$substance_use_any <- pmax(ahos1$tobacco_lifetime,
                               ahos1$tobacco_recent,
                               ahos1$alcohol_lifetime,
                               ahos1$alcohol_recent,
                               na.rm = TRUE)

ahos1$substance_use_any[is.infinite(ahos1$substance_use_any)] <- NA

table(ahos1$substance_use_any, useNA = "always")
prop.table(table(ahos1$substance_use_any, useNA = "always"))

#Alcohol/tobacco use by group#
table(ahos1$substance_use_any, ahos1$group, useNA = "always")
tab <- table(ahos1$substance_use_any, ahos1$group, useNA = "always")
prop.table(tab, margin = 2)
round(prop.table(tab, margin = 2) * 100, 1)

chisq.test(table(ahos1$substance_use_any, ahos1$group))
fisher.test(table(ahos1$substance_use_any, ahos1$group))

#Recent number of drinks#
table(ahos1$`_27`,useNA = "always")

#Lifetime number of drinks#
table(ahos1$`_28`,useNA = "always")

##Other drugs#
#Age of first drug usage#
ahos1$`_30`[ahos1$`_30` == 888] <- NA
table(ahos1$`_30`,useNA = "always")

#Lifetime bhang use#
ahos1$`_31`[ahos1$`_31` == 888] <- NA
table(ahos1$`_31`,useNA = "always")

#Lifetime miraa use#
ahos1$`_32`[ahos1$`_32` == 888] <- NA
table(ahos1$`_32`,useNA = "always")

#Recent bhang use#
ahos1$`_33`[ahos1$`_33` == 888] <- NA
table(ahos1$`_33`,useNA = "always")

#Recent miraa use#
ahos1$`_34`[ahos1$`_34` == 888] <- NA
table(ahos1$`_34`,useNA = "always")

#Lifetime glue use#
ahos1$`_35`[ahos1$`_35` == 888] <- NA
table(ahos1$`_35`,useNA = "always")



##########################################################################
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
Attention_items <- c("cbcl_sch_1","cbcl_sch_4","cbcl_sch_8","cbcl_sch_10","cbcl_sch_13",
                     "cbcl_sch_17","cbcl_sch_41","cbcl_sch_61","cbc_sch_l78","cbcl_sch_80")
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

by(
  data_child$cbcl_total_problems,
  data_child$group,
  summary
)

###Checking distribution#
# Histogram
hist(
  ahos1$cbcl_total_problems,
  breaks = 20,
  main = "CBCL Total Problems Score",
  xlab = "Total Problems Score"
)

# Density plot
plot(
  density(ahos1$cbcl_total_problems, na.rm = TRUE),
  main = "Density Plot of CBCL Total Problems Score"
)

# Summary statistics
summary(ahos1$cbcl_total_problems)

# Skewness and kurtosis
library(moments)

skewness(ahos1$cbcl_total_problems, na.rm = TRUE)
kurtosis(ahos1$cbcl_total_problems, na.rm = TRUE)


##Distribution by group#
library(ggplot2)

ggplot(ahos1, aes(x = cbcl_total_problems)) +
  geom_histogram(bins = 20) +
  facet_wrap(~group) +
  theme_bw()


ggplot(ahos1, aes(x = group, y = cbcl_total_problems)) +
  geom_boxplot() +
  theme_bw()

##Normality assessment#
# Overall
shapiro.test(ahos1$cbcl_total_problems)

# Within groups
by(
  ahos1$cbcl_total_problems,
  ahos1$group,
  shapiro.test
)

##Q-Q Plots#
qqnorm(ahos1$cbcl_total_problems)
qqline(ahos1$cbcl_total_problems)

#By group#
par(mfrow = c(1,3))

for(g in unique(ahos1$group)){
  x <- ahos1$cbcl_total_problems[ahos1$group == g]
  qqnorm(x, main = g)
  qqline(x)
}

par(mfrow = c(1,1))



####Now ANOVA and ANCOVA#
ahos1$group <- factor(
  ahos1$group,
  levels = c("control", "exposed", "infected")
)

anova_model <- aov(cbcl_total_problems ~ group, data = ahos1)
summary(anova_model)

## Bonferroni
#bonferroni = glht(anova_model, linfct = mcp(group = all_pairwise))#
#summary(bonferroni, test = adjusted("bonferroni"))#

library(multcomp)

bonferroni <- glht(
  anova_model,
  linfct = mcp(group = "Tukey")   # or "Bonferroni" via adjust argument
)

summary(bonferroni)

bonferroni_ci = confint(bonferroni)
bonferroni_ci
plot(bonferroni)

##ANCOVA for Total problems while adjusting for covariates##
#baseline_vars <- ahos[, c(
  "studyno",
  "depression",
  "sch_exp",
  "peer_exp",
  "medical_history",
  "disease_severity"
)]#

#ahos1_full <- merge(
  ahos1,
  baseline_vars,
  by = "id",
  all.x = TRUE
)#


##Main analysis looking at group differences#

library(car)

mod_total <- lm(
  cbcl_total_problems ~ group + csex + cage + child_edu + parents_alive1 + SES
  + sch_exp + peer_exp + medical_history + depression + stunted,
  data = ahos1
)


Anova(mod_total, type = 3)

summary(mod_total)

#Bonferroni adjusted pairwise comparisons#

library(multcomp)

summary(
  glht(mod_total, linfct = mcp(group = "Tukey")),
  test = adjusted("bonferroni")
)


###HIV-specific analysis#
##Examines effect of viral suppression vs non-suppression#
infected <- subset(ahos1, group == "infected")

infected$disease_severity <- relevel(
  factor(infected$disease_severity),
  ref = "Suppressed"
)

mod_inf <- lm(
  cbcl_total_problems ~ disease_severity + cage + child_edu +
    parents_alive1 + SES + sch_exp + peer_exp + csex +
    medical_history + depression + stunted,
  data = infected
)

Anova(mod_inf, type = 3)
summary(mod_inf )


##############################################
#Compare 4 groups#
#control
#exposed
#suppressed
#non-suppressed
############################################

#Create combined variable
library(dplyr)

ahos1 <- ahos1 %>%
  mutate(
    disease_group = case_when(
      
      group == "control" ~ "Control",
      group == "exposed" ~ "Exposed",
      
      group == "infected" & disease_severity == "Suppressed" ~ "Suppressed",
      group == "infected" & disease_severity == "Non-suppressed" ~ "Non-suppressed",
      
      TRUE ~ NA_character_
    )
  )

#Convert to factor#
ahos1$disease_group <- factor(
  ahos1$disease_group,
  levels = c("Control", "Exposed", "Suppressed", "Non-suppressed")
)

#Fit ANCOVA#
mod_ancova <- lm(
  cbcl_total_problems ~ disease_group + cage + child_edu +
    parents_alive1 + SES + sch_exp + peer_exp +
    medical_history + depression + stunted,
  data = ahos1
)

summary(mod_ancova)


##Model diagnostics#
# Normality of residuals
qqnorm(residuals(mod_total))
qqline(residuals(mod_total))

# Homogeneity of variance
leveneTest(TotalProbs  ~ group2, data = ahos1)


##MANOVA for subscales##
Anxious_items <- c("cbcl14","cbcl29","cbcl30","cbcl31",
                   "cbcl32","cbcl33","cbcl35","cbcl45",
                   "cbcl50","cbcl52","cbcl71","cbcl91","cbcl112")
Withdrawn_items <- c("cbcl5","cbcl42","cbcl65","cbcl69",
                     "cbcl75","cbcl102","cbcl103","cbcl111")
Somatic_items <- c("cbcl47","cbcl49","cbcl51","cbcl54",
                   "cbcl56a","cbcl56b","cbcl56c","cbcl56d",
                   "cbcl56e","cbcl56f","cbcl56g")
Social_items <- c("cbcl11","cbcl12","cbcl25","cbcl27",
                  "cbcl34","cbcl36","cbcl38","cbcl48",
                  "cbcl62","cbcl64","cbcl79")
Thought_items <- c("cbcl9","cbcl18","cbcl40","cbcl46","cbcl58","cbcl59","cbcl60",
                   "cbcl66","cbcl70","cbcl76","cbcl83","cbcl84","cbcl85","cbcl92","cbcl100")
Attention_items <- c("cbcl1","cbcl4","cbcl8","cbcl10","cbcl13",
                     "cbcl17","cbcl41","cbcl61","cbcl78","cbcl80")
Rule_items <- c("cbcl2","cbcl26","cbcl28","cbcl39","cbcl43",
                "cbcl63","cbcl67","cbcl72","cbcl73","cbcl81",
                "cbcl82","cbcl90","cbcl96","cbcl99",
                "cbcl101","cbcl105","cbcl106")
Aggressive_items <- c("cbcl3","cbcl6","cbcl19","cbcl20",
                      "cbcl21","cbcl22","cbcl23","cbcl37",
                      "cbcl57","cbcl68","cbcl86","cbcl87",
                      "cbcl88","cbcl89","cbcl94","cbcl95",
                      "cbcl97","cbcl104")
Other_items <- c("cbcl6","cbcl7","cbcl15","cbcl24","cbcl44",
                 "cbcl53","cbcl55","cbcl74","cbcl77","cbcl93",
                 "cbcl98","cbcl107","cbcl108","cbcl109","cbcl110")
Internalizing_items <- c("cbcl14","cbcl29","cbcl30","cbcl31",
                         "cbcl32","cbcl33","cbcl35","cbcl45",
                         "cbcl50","cbcl52","cbcl71","cbcl91","cbcl112",
                         "cbcl5","cbcl42","cbcl65","cbcl69",
                         "cbcl75","cbcl102","cbcl103","cbcl111",
                         "cbcl47","cbcl49","cbcl51","cbcl54",
                         "cbcl56a","cbcl56b","cbcl56c","cbcl56d",
                         "cbcl56e","cbcl56f","cbcl56g")
Externalizing_items <- c("cbcl2","cbcl26","cbcl28","cbcl39","cbcl43",
                         "cbcl63","cbcl67","cbcl72","cbcl81",
                         "cbcl82","cbcl90","cbcl96","cbcl99",
                         "cbcl101","cbcl105","cbcl106",
                         "cbcl3","cbcl6","cbcl19","cbcl20",
                         "cbcl21","cbcl22","cbcl23","cbcl37",
                         "cbcl57","cbcl68","cbcl86","cbcl87",
                         "cbcl88","cbcl89","cbcl94","cbcl95",
                         "cbcl97","cbcl104")

##Compute row totals##
ahos1 <- ahos1 %>%
  mutate(
    anxious_depressed = rowSums(across(all_of(Anxious_items)), na.rm = TRUE),
    withdrawn_depressed = rowSums(across(all_of(Withdrawn_items)), na.rm = TRUE),
    somatic_complaints = rowSums(across(all_of(Somatic_items)), na.rm = TRUE),
    social_problems = rowSums(across(all_of(Social_items)), na.rm = TRUE),
    thought_problems = rowSums(across(all_of(Thought_items)), na.rm = TRUE),
    attention_problems = rowSums(across(all_of(Attention_items)), na.rm = TRUE),
    rule_breaking = rowSums(across(all_of(Rule_items)), na.rm = TRUE),
    aggressive_behavior = rowSums(across(all_of(Aggressive_items)), na.rm = TRUE),
    other_problems = rowSums(across(all_of(Other_items)), na.rm = TRUE),
    internalizing_problems = rowSums(across(all_of(Internalizing_items)), na.rm = TRUE),
    externalizing_problems = rowSums(across(all_of(Externalizing_items)), na.rm = TRUE)
  )

cbcl_subscales <- ahos1 %>%
  select(
    group,
    child_edu,
    cage,
    csex,
    parents_alive1,
    SES,
    sch_exp,
    peer_exp,
    medical_history,
    depression,
    stunted,
    anxious_depressed,
    withdrawn_depressed,
    somatic_complaints,
    social_problems,
    thought_problems,
    attention_problems,
    rule_breaking,
    aggressive_behavior,
    other_problems,
    internalizing_problems,
    externalizing_problems
  ) %>%
  drop_na()

#Fit manova#
manova_mod <- manova(
  cbind(
    anxious_depressed,
    withdrawn_depressed,
    somatic_complaints,
    social_problems,
    thought_problems,
    attention_problems,
    rule_breaking,
    aggressive_behavior
  ) ~ group,
  data = cbcl_subscales
)

#Test results#
summary(manova_mod, test = "Pillai")

summary.aov(manova_mod)

#Bonferroni correction#
library(multcomp)

subscales <- c(
  "anxious_depressed",
  "withdrawn_depressed",
  "somatic_complaints",
  "social_problems",
  "thought_problems",
  "attention_problems",
  "rule_breaking",
  "aggressive_behavior"
)

bonf_results <- lapply(subscales, function(s) {
  mod <- aov(as.formula(paste(s, "~ group")), data = cbcl_subscales)
  summary(
    glht(mod, linfct = mcp(group = "Tukey")),
    test = adjusted("bonferroni")
  )
})

names(bonf_results) <- subscales

##Extracting bonferroni results#
bonf_results$attention_problems
bonf_results$anxious_depressed
bonf_results$withdrawn_depressed
bonf_results$somatic_complaints
bonf_results$social_problems
bonf_results$thought_problems
bonf_results$rule_breaking
bonf_results$aggressive_behavior

#MANOVA assumptions#

library(biotools)

boxM(                    #homogeneity of covariance matrices
  cbcl_subscales[, -1],  # remove group column
  cbcl_subscales$group2
)

cor(cbcl_subscales[, -1], use = "pairwise.complete.obs")


##MANOVA for broadband scales##

cbcl_subscales1 <- ahos1 %>%
  select(
    group,
    internalizing_problems,
    externalizing_problems
  ) %>%
  drop_na()

#Fit manova#
manova_mod1 <- manova(
  cbind(
    internalizing_problems,
    externalizing_problems
  ) ~ group,
  data = cbcl_subscales1
)

#Test results#
summary(manova_mod1, test = "Pillai")

summary.aov(manova_mod1)

#Bonferroni correction#
library(multcomp)

subscales1 <- c(
  "internalizing_problems",
  "externalizing_problems"
)

bonf_results1 <- lapply(subscales1, function(s) {
  mod <- aov(as.formula(paste(s, "~ group")), data = cbcl_subscales1)
  summary(
    glht(mod, linfct = mcp(group = "Tukey")),
    test = adjusted("bonferroni")
  )
})

names(bonf_results1) <- subscales1

##Extracting bonferroni results#
bonf_results1$internalizing_problems
bonf_results1$externalizing_problems


#MANOVA assumptions#

library(biotools)

boxM(                    #homogeneity of covariance matrices
  cbcl_subscales[, -1],  # remove group column
  cbcl_subscales$group2
)

cor(cbcl_subscales[, -1], use = "pairwise.complete.obs")

##MANCOVA (adding covariates)#
manova_cov <- manova(
  cbind(
    anxious_depressed,
    withdrawn_depressed,
    somatic_complaints,
    social_problems,
    thought_problems,
    attention_problems,
    rule_breaking,
    aggressive_behavior
  ) ~ group + cage + csex +child_edu + parents_alive1 + SES + sch_exp + 
    peer_exp + medical_history + depression + stunted,
  data = cbcl_subscales
)

summary(manova_cov, test = "Pillai")
summary.aov(manova_cov)

##Covariate-adjusted univariate ANCOVAs#
library(car)

subscales <- c(
  "anxious_depressed",
  "withdrawn_depressed",
  "somatic_complaints",
  "social_problems",
  "thought_problems",
  "attention_problems",
  "rule_breaking",
  "aggressive_behavior"
)

ancova_results <- lapply(subscales, function(s) {
  mod <- lm(
    as.formula(
      paste(
        s,
        "~ group2 + csex + SES + depression"
      )
    ),
    data = cbcl_subscales
  )
  Anova(mod, type = 3)
})

names(ancova_results) <- subscales

##Extracting covariate-adjusted ANCOVA results#
ancova_results$attention_problems
ancova_results$anxious_depressed
ancova_results$withdrawn_depressed
ancova_results$somatic_complaints
ancova_results$social_problems
ancova_results$thought_problems
ancova_results$rule_breaking
ancova_results$aggressive_behavior

##Bonferroni-adjusted pairwise comparisons (HIV status)#
library(multcomp)

pairwise_results <- lapply(subscales, function(s) {
  mod <- lm(
    as.formula(
      paste(
        s,
        "~ group2 + csex + SES + depression"
      )
    ),
    data = cbcl_subscales
  )
  
  summary(
    glht(mod, linfct = mcp(group2 = "Tukey")),
    test = adjusted("bonferroni")
  )
})

names(pairwise_results) <- subscales

##Extracting bonferroni results#
pairwise_results$attention_problems
pairwise_results$anxious_depressed
pairwise_results$withdrawn_depressed
pairwise_results$somatic_complaints
pairwise_results$social_problems
pairwise_results$thought_problems
pairwise_results$rule_breaking
pairwise_results$aggressive_behavior

##Export tables to word##
library(broom)
library(dplyr)
library(purrr)
library(flextable)
library(officer)
library(multcomp)
library(car)

subscales <- c(
  "anxious_depressed",
  "withdrawn_depressed",
  "somatic_complaints",
  "social_problems",
  "thought_problems",
  "attention_problems",
  "rule_breaking",
  "aggressive_behavior"
)

ancova_tbl <- map_dfr(subscales, function(s) {
  mod <- lm(
    as.formula(
      paste(s, "~ group2 + child_edu + parents_alive1 + SES + depression")
    ),
    data = cbcl_subscales
  )
  
  tidy(Anova(mod, type = 3)) %>%
    filter(term == "group2") %>%
    mutate(Subscale = s) %>%
    select(Subscale, df, statistic, p.value)
})

ancova_tbl <- ancova_tbl %>%
  mutate(
    p.value = ifelse(p.value < .001, "< .001", round(p.value, 3)),
    F = round(statistic, 2)
  ) %>%
  select(Subscale, df, F, p.value)

ancova_tbl

bonf_tbl <- map_dfr(subscales, function(s) {
  mod <- lm(
    as.formula(
      paste(s, "~ group2 + child_edu + parents_alive1 + SES + depression")
    ),
    data = cbcl_subscales
  )
  
  glht_res <- summary(
    glht(mod, linfct = mcp(group2 = "Tukey")),
    test = adjusted("bonferroni")
  )
  
  as.data.frame(glht_res$test[c("coefficients", "pvalues")]) %>%
    mutate(
      Contrast = rownames(.),
      Subscale = s,
      Mean_Diff = round(coefficients, 2),
      p.value = ifelse(pvalues < .001, "< .001", round(pvalues, 3))
    ) %>%
    select(Subscale, Contrast, Mean_Diff, p.value)
})



##MANCOVA for broadband scales#
manova_cov1 <- manova(
  cbind(
    internalizing_problems,
    externalizing_problems
  ) ~ group2 + cage + child_edu + parents_alive1 + SES + sch_exp + 
    peer_exp + medical_history + depression,
  data = cbcl_subscales
)

summary(manova_cov1, test = "Pillai")


##Covariate-adjusted univariate ANCOVAs#
library(car)

subscales <- c(
  "internalizing_problems",
  "externalizing_problems"
)

ancova_results1 <- lapply(subscales, function(s) {
  mod <- lm(
    as.formula(
      paste(
        s,
        "~ group2 + child_edu + parents_alive1 + SES + depression"
      )
    ),
    data = cbcl_subscales
  )
  Anova(mod, type = 3)
})

names(ancova_results1) <- subscales

##Extracting covariate-adjusted ANCOVA results#
ancova_results1$internalizing_problems
ancova_results1$externalizing_problems


##Bonferroni-adjusted pairwise comparisons (HIV status)#
library(multcomp)

pairwise_results1 <- lapply(subscales, function(s) {
  mod <- lm(
    as.formula(
      paste(
        s,
        "~ group2 + child_edu + parents_alive1 + SES + depression"
      )
    ),
    data = cbcl_subscales
  )
  
  summary(
    glht(mod, linfct = mcp(group2 = "Tukey")),
    test = adjusted("bonferroni")
  )
})

names(pairwise_results1) <- subscales

##Extracting bonferroni results#
pairwise_results1$internalizing_problems
pairwise_results1$externalizing_problems





####################################
#PAUL#
#####################################
library(multcomp)
ahos1$group2 = as.factor(ifelse(ahos1$group == "control",1,ifelse(ahos1$group == "exposed",2,3)))
table(ahos1$group2)
#mod <- aov(TotalProbs  ~group2,data = ahos1)
#tukey = glht(mod, linfct = mcp(group2 = "Tukey"))
#tukey_ci = confint(tukey)
#tukey_ci

#plot(tukey)

all_pairwise = rbind("2 - 1" = c(-1, 1, 0),
                     "3 - 1" = c(-1, 0, 1),
                     "3 - 2" = c(0, -1, 1))

mod <- aov(TotalProbs  ~group2,data = ahos1)
summary(mod)

## Bonferroni
bonferroni = glht(mod, linfct = mcp(group2 = all_pairwise))
summary(bonferroni, test = adjusted("bonferroni"))

bonferroni_ci = confint(bonferroni)
bonferroni_ci
plot(bonferroni)

##ANCOVA for Total problems while adjusting for covariates##
library(car)

mod_total <- lm(
  TotalProbs ~ group2 + cage + child_edu + parents_alive1 + SES + sch_exp + peer_exp + medical_history +depression,
  data = ahos1
)

mod_total <- lm(
  TotalProbs ~ group2 + cage + csex + parents_alive1 + SES + sch_exp + medical_history +depression,
  data = ahos1
)

Anova(mod_total, type = 3)

#Bonferroni adjusted pairwise comparisons#
library(multcomp)

summary(
  glht(mod_total, linfct = mcp(group2 = "Tukey")),
  test = adjusted("bonferroni")
)

##Model diagnostics#
# Normality of residuals
qqnorm(residuals(mod_total))
qqline(residuals(mod_total))

# Homogeneity of variance
leveneTest(TotalProbs  ~ group2, data = ahos1)


##MANOVA for subscales##
Anxious_items <- c("cbcl14","cbcl29","cbcl30","cbcl31",
                   "cbcl32","cbcl33","cbcl35","cbcl45",
                   "cbcl50","cbcl52","cbcl71","cbcl91","cbcl112")
Withdrawn_items <- c("cbcl5","cbcl42","cbcl65","cbcl69",
                     "cbcl75","cbcl102","cbcl103","cbcl111")
Somatic_items <- c("cbcl47","cbcl49","cbcl51","cbcl54",
                   "cbcl56a","cbcl56b","cbcl56c","cbcl56d",
                   "cbcl56e","cbcl56f","cbcl56g")
Social_items <- c("cbcl11","cbcl12","cbcl25","cbcl27",
                  "cbcl34","cbcl36","cbcl38","cbcl48",
                  "cbcl62","cbcl64","cbcl79")
Thought_items <- c("cbcl9","cbcl18","cbcl40","cbcl46","cbcl58","cbcl59","cbcl60",
                   "cbcl66","cbcl70","cbcl76","cbcl83","cbcl84","cbcl85","cbcl92","cbcl100")
Attention_items <- c("cbcl1","cbcl4","cbcl8","cbcl10","cbcl13",
                     "cbcl17","cbcl41","cbcl61","cbcl78","cbcl80")
Rule_items <- c("cbcl2","cbcl26","cbcl28","cbcl39","cbcl43",
                "cbcl63","cbcl67","cbcl72","cbcl73","cbcl81",
                "cbcl82","cbcl90","cbcl96","cbcl99",
                "cbcl101","cbcl105","cbcl106")
Aggressive_items <- c("cbcl3","cbcl6","cbcl19","cbcl20",
                      "cbcl21","cbcl22","cbcl23","cbcl37",
                      "cbcl57","cbcl68","cbcl86","cbcl87",
                      "cbcl88","cbcl89","cbcl94","cbcl95",
                      "cbcl97","cbcl104")
Other_items <- c("cbcl6","cbcl7","cbcl15","cbcl24","cbcl44",
                 "cbcl53","cbcl55","cbcl74","cbcl77","cbcl93",
                 "cbcl98","cbcl107","cbcl108","cbcl109","cbcl110")
Internalizing_items <- c("cbcl14","cbcl29","cbcl30","cbcl31",
                         "cbcl32","cbcl33","cbcl35","cbcl45",
                         "cbcl50","cbcl52","cbcl71","cbcl91","cbcl112",
                         "cbcl5","cbcl42","cbcl65","cbcl69",
                         "cbcl75","cbcl102","cbcl103","cbcl111",
                         "cbcl47","cbcl49","cbcl51","cbcl54",
                         "cbcl56a","cbcl56b","cbcl56c","cbcl56d",
                         "cbcl56e","cbcl56f","cbcl56g")
Externalizing_items <- c("cbcl2","cbcl26","cbcl28","cbcl39","cbcl43",
                         "cbcl63","cbcl67","cbcl72","cbcl81",
                         "cbcl82","cbcl90","cbcl96","cbcl99",
                         "cbcl101","cbcl105","cbcl106",
                         "cbcl3","cbcl6","cbcl19","cbcl20",
                         "cbcl21","cbcl22","cbcl23","cbcl37",
                         "cbcl57","cbcl68","cbcl86","cbcl87",
                         "cbcl88","cbcl89","cbcl94","cbcl95",
                         "cbcl97","cbcl104")

##Compute row totals##
ahos1 <- ahos1 %>%
  mutate(
    anxious_depressed = rowSums(across(all_of(Anxious_items)), na.rm = TRUE),
    withdrawn_depressed = rowSums(across(all_of(Withdrawn_items)), na.rm = TRUE),
    somatic_complaints = rowSums(across(all_of(Somatic_items)), na.rm = TRUE),
    social_problems = rowSums(across(all_of(Social_items)), na.rm = TRUE),
    thought_problems = rowSums(across(all_of(Thought_items)), na.rm = TRUE),
    attention_problems = rowSums(across(all_of(Attention_items)), na.rm = TRUE),
    rule_breaking = rowSums(across(all_of(Rule_items)), na.rm = TRUE),
    aggressive_behavior = rowSums(across(all_of(Aggressive_items)), na.rm = TRUE),
    other_problems = rowSums(across(all_of(Other_items)), na.rm = TRUE),
    internalizing_problems = rowSums(across(all_of(Internalizing_items)), na.rm = TRUE),
    externalizing_problems = rowSums(across(all_of(Externalizing_items)), na.rm = TRUE)
  )

cbcl_subscales <- ahos1 %>%
  select(
    group2,
    child_edu,
    cage,
    csex,
    parents_alive1,
    SES,
    sch_exp,
    peer_exp,
    medical_history,
    depression,
    anxious_depressed,
    withdrawn_depressed,
    somatic_complaints,
    social_problems,
    thought_problems,
    attention_problems,
    rule_breaking,
    aggressive_behavior,
    other_problems,
    internalizing_problems,
    externalizing_problems
  ) %>%
  drop_na()

#Fit manova#
manova_mod <- manova(
  cbind(
    anxious_depressed,
    withdrawn_depressed,
    somatic_complaints,
    social_problems,
    thought_problems,
    attention_problems,
    rule_breaking,
    aggressive_behavior
  ) ~ group2,
  data = cbcl_subscales
)

#Test results#
summary(manova_mod, test = "Pillai")

summary.aov(manova_mod)

#Bonferroni correction#
library(multcomp)

subscales <- c(
  "anxious_depressed",
  "withdrawn_depressed",
  "somatic_complaints",
  "social_problems",
  "thought_problems",
  "attention_problems",
  "rule_breaking",
  "aggressive_behavior"
)

bonf_results <- lapply(subscales, function(s) {
  mod <- aov(as.formula(paste(s, "~ group2")), data = cbcl_subscales)
  summary(
    glht(mod, linfct = mcp(group2 = "Tukey")),
    test = adjusted("bonferroni")
  )
})

names(bonf_results) <- subscales

##Extracting bonferroni results#
bonf_results$attention_problems
bonf_results$anxious_depressed
bonf_results$withdrawn_depressed
bonf_results$somatic_complaints
bonf_results$social_problems
bonf_results$thought_problems
bonf_results$rule_breaking
bonf_results$aggressive_behavior

#MANOVA assumptions#

library(biotools)

boxM(                    #homogeneity of covariance matrices
  cbcl_subscales[, -1],  # remove group column
  cbcl_subscales$group2
)

cor(cbcl_subscales[, -1], use = "pairwise.complete.obs")


##MANOVA for broadband scales##

cbcl_subscales1 <- ahos1 %>%
  select(
    group2,
    internalizing_problems,
    externalizing_problems
  ) %>%
  drop_na()

#Fit manova#
manova_mod1 <- manova(
  cbind(
    internalizing_problems,
    externalizing_problems
  ) ~ group2,
  data = cbcl_subscales1
)

#Test results#
summary(manova_mod1, test = "Pillai")

summary.aov(manova_mod1)

#Bonferroni correction#
library(multcomp)

subscales1 <- c(
  "internalizing_problems",
  "externalizing_problems"
)

bonf_results1 <- lapply(subscales1, function(s) {
  mod <- aov(as.formula(paste(s, "~ group2")), data = cbcl_subscales1)
  summary(
    glht(mod, linfct = mcp(group2 = "Tukey")),
    test = adjusted("bonferroni")
  )
})

names(bonf_results1) <- subscales1

##Extracting bonferroni results#
bonf_results1$internalizing_problems
bonf_results1$externalizing_problems


#MANOVA assumptions#

library(biotools)

boxM(                    #homogeneity of covariance matrices
  cbcl_subscales[, -1],  # remove group column
  cbcl_subscales$group2
)

cor(cbcl_subscales[, -1], use = "pairwise.complete.obs")

##MANCOVA (adding covariates)#
manova_cov <- manova(
  cbind(
    anxious_depressed,
    withdrawn_depressed,
    somatic_complaints,
    social_problems,
    thought_problems,
    attention_problems,
    rule_breaking,
    aggressive_behavior
  ) ~ group2 + cage + child_edu + parents_alive1 + SES + sch_exp + 
    peer_exp + medical_history + depression,
  data = cbcl_subscales
)

manova_cov <- manova(
  cbind(
    anxious_depressed,
    withdrawn_depressed,
    somatic_complaints,
    social_problems,
    thought_problems,
    attention_problems,
    rule_breaking,
    aggressive_behavior
  ) ~ group2 + cage + csex + parents_alive1 + SES + sch_exp + 
    medical_history + depression,
  data = cbcl_subscales
)

summary(manova_cov, test = "Pillai")


##Covariate-adjusted univariate ANCOVAs#
library(car)

subscales <- c(
  "anxious_depressed",
  "withdrawn_depressed",
  "somatic_complaints",
  "social_problems",
  "thought_problems",
  "attention_problems",
  "rule_breaking",
  "aggressive_behavior"
)

ancova_results <- lapply(subscales, function(s) {
  mod <- lm(
    as.formula(
      paste(
        s,
        "~ group2 + csex + SES + depression"
      )
    ),
    data = cbcl_subscales
  )
  Anova(mod, type = 3)
})

names(ancova_results) <- subscales

##Extracting covariate-adjusted ANCOVA results#
ancova_results$attention_problems
ancova_results$anxious_depressed
ancova_results$withdrawn_depressed
ancova_results$somatic_complaints
ancova_results$social_problems
ancova_results$thought_problems
ancova_results$rule_breaking
ancova_results$aggressive_behavior

##Bonferroni-adjusted pairwise comparisons (HIV status)#
library(multcomp)

pairwise_results <- lapply(subscales, function(s) {
  mod <- lm(
    as.formula(
      paste(
        s,
        "~ group2 + csex + SES + depression"
      )
    ),
    data = cbcl_subscales
  )
  
  summary(
    glht(mod, linfct = mcp(group2 = "Tukey")),
    test = adjusted("bonferroni")
  )
})

names(pairwise_results) <- subscales

##Extracting bonferroni results#
pairwise_results$attention_problems
pairwise_results$anxious_depressed
pairwise_results$withdrawn_depressed
pairwise_results$somatic_complaints
pairwise_results$social_problems
pairwise_results$thought_problems
pairwise_results$rule_breaking
pairwise_results$aggressive_behavior

##Export tables to word##
library(broom)
library(dplyr)
library(purrr)
library(flextable)
library(officer)
library(multcomp)
library(car)

subscales <- c(
  "anxious_depressed",
  "withdrawn_depressed",
  "somatic_complaints",
  "social_problems",
  "thought_problems",
  "attention_problems",
  "rule_breaking",
  "aggressive_behavior"
)

ancova_tbl <- map_dfr(subscales, function(s) {
  mod <- lm(
    as.formula(
      paste(s, "~ group2 + child_edu + parents_alive1 + SES + depression")
    ),
    data = cbcl_subscales
  )
  
  tidy(Anova(mod, type = 3)) %>%
    filter(term == "group2") %>%
    mutate(Subscale = s) %>%
    select(Subscale, df, statistic, p.value)
})

ancova_tbl <- ancova_tbl %>%
  mutate(
    p.value = ifelse(p.value < .001, "< .001", round(p.value, 3)),
    F = round(statistic, 2)
  ) %>%
  select(Subscale, df, F, p.value)

ancova_tbl

bonf_tbl <- map_dfr(subscales, function(s) {
  mod <- lm(
    as.formula(
      paste(s, "~ group2 + child_edu + parents_alive1 + SES + depression")
    ),
    data = cbcl_subscales
  )
  
  glht_res <- summary(
    glht(mod, linfct = mcp(group2 = "Tukey")),
    test = adjusted("bonferroni")
  )
  
  as.data.frame(glht_res$test[c("coefficients", "pvalues")]) %>%
    mutate(
      Contrast = rownames(.),
      Subscale = s,
      Mean_Diff = round(coefficients, 2),
      p.value = ifelse(pvalues < .001, "< .001", round(pvalues, 3))
    ) %>%
    select(Subscale, Contrast, Mean_Diff, p.value)
})



##MANCOVA for broadband scales#
manova_cov1 <- manova(
  cbind(
    internalizing_problems,
    externalizing_problems
  ) ~ group2 + cage + child_edu + parents_alive1 + SES + sch_exp + 
    peer_exp + medical_history + depression,
  data = cbcl_subscales
)

summary(manova_cov1, test = "Pillai")


##Covariate-adjusted univariate ANCOVAs#
library(car)

subscales <- c(
  "internalizing_problems",
  "externalizing_problems"
)

ancova_results1 <- lapply(subscales, function(s) {
  mod <- lm(
    as.formula(
      paste(
        s,
        "~ group2 + child_edu + parents_alive1 + SES + depression"
      )
    ),
    data = cbcl_subscales
  )
  Anova(mod, type = 3)
})

names(ancova_results1) <- subscales

##Extracting covariate-adjusted ANCOVA results#
ancova_results1$internalizing_problems
ancova_results1$externalizing_problems


##Bonferroni-adjusted pairwise comparisons (HIV status)#
library(multcomp)

pairwise_results1 <- lapply(subscales, function(s) {
  mod <- lm(
    as.formula(
      paste(
        s,
        "~ group2 + child_edu + parents_alive1 + SES + depression"
      )
    ),
    data = cbcl_subscales
  )
  
  summary(
    glht(mod, linfct = mcp(group2 = "Tukey")),
    test = adjusted("bonferroni")
  )
})

names(pairwise_results1) <- subscales

##Extracting bonferroni results#
pairwise_results1$internalizing_problems
pairwise_results1$externalizing_problems


#Mean scores and bar graph##
# Load libraries
library(tidyverse)

# Example: load your dataset
cbcl_all <- read.csv("Final cbcl analysis data.csv")

# ---- Step 1. Define group variable ----
# e.g., group variable could be "status" with values "Infected", "Exposed", "Community control"
group_var <- "group"  


# Example: define your CBCL items per subscale
Anxious_items <- c("cbcl14","cbcl29","cbcl30","cbcl31",
                   "cbcl32","cbcl33","cbcl35","cbcl45",
                   "cbcl50","cbcl52","cbcl71","cbcl91","cbcl112")
Withdrawn_items <- c("cbcl5","cbcl42","cbcl65","cbcl69",
                     "cbcl75","cbcl102","cbcl103","cbcl111")
Somatic_items <- c("cbcl47","cbcl49","cbcl51","cbcl54",
                   "cbcl56a","cbcl56b","cbcl56c","cbcl56d",
                   "cbcl56e","cbcl56f","cbcl56g")
Social_items <- c("cbcl11","cbcl12","cbcl25","cbcl27",
                  "cbcl34","cbcl36","cbcl38","cbcl48",
                  "cbcl62","cbcl64","cbcl79")
Thought_items <- c("cbcl9","cbcl18","cbcl40","cbcl46","cbcl58","cbcl59","cbcl60",
                   "cbcl66","cbcl70","cbcl76","cbcl83","cbcl84","cbcl85","cbcl92","cbcl100")
Attention_items <- c("cbcl1","cbcl4","cbcl8","cbcl10","cbcl13",
                     "cbcl17","cbcl41","cbcl61","cbcl78","cbcl80")
Rule_items <- c("cbcl2","cbcl26","cbcl28","cbcl39","cbcl43",
                "cbcl63","cbcl67","cbcl72","cbcl81",
                "cbcl82","cbcl90","cbcl96","cbcl99",
                "cbcl101","cbcl105","cbcl106")
Aggressive_items <- c("cbcl3","cbcl6","cbcl19","cbcl20",
                      "cbcl21","cbcl22","cbcl23","cbcl37",
                      "cbcl57","cbcl68","cbcl86","cbcl87",
                      "cbcl88","cbcl89","cbcl94","cbcl95",
                      "cbcl97","cbcl104")
Other_items <- c("cbcl6","cbcl7","cbcl15","cbcl24","cbcl44",
                 "cbcl53","cbcl55","cbcl74","cbcl77","cbcl93",
                 "cbcl98","cbcl107","cbcl108","cbcl109","cbcl110")
Internalizing_items <- c("cbcl14","cbcl29","cbcl30","cbcl31",
                         "cbcl32","cbcl33","cbcl35","cbcl45",
                         "cbcl50","cbcl52","cbcl71","cbcl91","cbcl112",
                         "cbcl5","cbcl42","cbcl65","cbcl69",
                         "cbcl75","cbcl102","cbcl103","cbcl111",
                         "cbcl47","cbcl49","cbcl51","cbcl54",
                         "cbcl56a","cbcl56b","cbcl56c","cbcl56d",
                         "cbcl56e","cbcl56f","cbcl56g")
Externalizing_items <- c("cbcl2","cbcl26","cbcl28","cbcl39","cbcl43",
                         "cbcl63","cbcl67","cbcl72","cbcl81",
                         "cbcl82","cbcl90","cbcl96","cbcl99",
                         "cbcl101","cbcl105","cbcl106",
                         "cbcl3","cbcl6","cbcl19","cbcl20",
                         "cbcl21","cbcl22","cbcl23","cbcl37",
                         "cbcl57","cbcl68","cbcl86","cbcl87",
                         "cbcl88","cbcl89","cbcl94","cbcl95",
                         "cbcl97","cbcl104")

##Compute row totals##
cbcl_all <- cbcl_all %>%
  mutate(
    anxious_depressed = rowSums(across(all_of(Anxious_items)), na.rm = TRUE),
    withdrawn_depressed = rowSums(across(all_of(Withdrawn_items)), na.rm = TRUE),
    somatic_complaints = rowSums(across(all_of(Somatic_items)), na.rm = TRUE),
    social_problems = rowSums(across(all_of(Social_items)), na.rm = TRUE),
    thought_problems = rowSums(across(all_of(Thought_items)), na.rm = TRUE),
    attention_problems = rowSums(across(all_of(Attention_items)), na.rm = TRUE),
    rule_breaking = rowSums(across(all_of(Rule_items)), na.rm = TRUE),
    aggressive_behavior = rowSums(across(all_of(Aggressive_items)), na.rm = TRUE),
    other_problems = rowSums(across(all_of(Other_items)), na.rm = TRUE),
    internalizing_problems = rowSums(across(all_of(Internalizing_items)), na.rm = TRUE),
    externalizing_problems = rowSums(across(all_of(Externalizing_items)), na.rm = TRUE)
  )

##Compute means##
cbcl_vars <- c("anxious_depressed", "withdrawn_depressed", "somatic_complaints",
               "social_problems", "thought_problems", "attention_problems",
               "rule_breaking", "aggressive_behavior", "other_problems",
               "internalizing_problems", "externalizing_problems")  

cbcl_means_total <- cbcl_all %>%
  summarise(across(all_of(cbcl_vars),
                   list(mean = ~mean(.x, na.rm = TRUE),
                        sd = ~sd(.x, na.rm = TRUE)),
                   .names = "{.col}_{.fn}"))

cbcl_means_total


cbcl_means_by_group <- cbcl_all %>%
  group_by(group) %>%
  summarise(across(all_of(cbcl_vars), ~ mean(.x, na.rm = TRUE)))

cbcl_means_by_group

cbcl_means_by_sex <- cbcl_all %>%
  group_by(csex) %>%
  summarise(across(all_of(cbcl_vars), ~ mean(.x, na.rm = TRUE)))

cbcl_means_by_sex



##CBCL subscale scoring##
library(dplyr)

# Function to compute mean-of-means, SD-of-means, SE-of-means
subscale_stats <- function(df, items) {
  
  # Compute item means, sds, ses
  stats <- df %>%
    summarise(across(
      all_of(items),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE),
        se   = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
      ),
      .names = "{.col}_{.fn}"
    ))
  
  # Extract item means
  item_means <- stats %>% select(ends_with("_mean"))
  
  # Compute subscale mean-of-means
  subscale_mean <- mean(as.numeric(item_means))
  
  # SD of item means
  subscale_sd <- sd(as.numeric(item_means))
  
  # SE = SD / sqrt(n_items)
  subscale_se <- subscale_sd / sqrt(ncol(item_means))
  
  return(data.frame(
    mean = subscale_mean,
    sd = subscale_sd,
    se = subscale_se
  ))
}

##Whole-sample subscale scores##
cbcl_whole <- bind_rows(
  Withdrawn      = subscale_stats(cbcl_all, Withdrawn_items),
  Anxious        = subscale_stats(cbcl_all, Anxious_items),
  Somatic        = subscale_stats(cbcl_all, Somatic_items),
  Social         = subscale_stats(cbcl_all, Social_items),
  Thought        = subscale_stats(cbcl_all, Thought_items),
  Attention      = subscale_stats(cbcl_all, Attention_items),
  RuleBreaking   = subscale_stats(cbcl_all, Rule_items),
  Aggressive     = subscale_stats(cbcl_all, Aggressive_items),
  Other          = subscale_stats(cbcl_all, Other_items),
  Internalizing  = subscale_stats(cbcl_all, Internalizing_items),
  Externalizing  = subscale_stats(cbcl_all, Externalizing_items),
  .id = "subscale"
)

cbcl_whole

##Subscale scores by group##
groups <- unique(cbcl_all$group)

cbcl_by_group <- lapply(groups, function(g) {
  
  df_group <- cbcl_all %>% filter(group == g)
  
  out <- bind_rows(
    Withdrawn      = subscale_stats(df_group, Withdrawn_items),
    Anxious        = subscale_stats(df_group, Anxious_items),
    Somatic        = subscale_stats(df_group, Somatic_items),
    Social         = subscale_stats(df_group, Social_items),
    Thought        = subscale_stats(df_group, Thought_items),
    Attention      = subscale_stats(df_group, Attention_items),
    RuleBreaking   = subscale_stats(df_group, Rule_items),
    Aggressive     = subscale_stats(df_group, Aggressive_items),
    Other          = subscale_stats(df_group, Other_items),
    Internalizing  = subscale_stats(df_group, Internalizing_items),
    Externalizing  = subscale_stats(df_group, Externalizing_items),
    .id = "subscale"
  )
  
  out$group <- g
  return(out)
})

cbcl_by_group <- bind_rows(cbcl_by_group)

cbcl_by_group


##Confirming mean,Doro##

cbcl_items=cbcl_all[, c("cbcl14", "cbcl29", "cbcl30", "cbcl31", "cbcl32",
                        "cbcl33", "cbcl35", "cbcl45", "cbcl50", "cbcl52",
                        "cbcl71", "cbcl91", "cbcl112")]

item_means=colMeans(cbcl_items, na.rm = T)
print(item_means)
sum(item_means)

##Plot with derived means##
library(ggplot2)
library(dplyr)

# Ensure subscales appear in a nice consistent order
cbcl_by_group$subscale <- factor(cbcl_by_group$subscale,
                                 levels = unique(cbcl_by_group$subscale))

# Plot
ggplot(cbcl_by_group, aes(x = subscale, y = mean, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.8),
                width = 0.2) +
  coord_flip() +
  labs(
    x = "CBCL Subscale",
    y = "Mean Score (mean of item means)",
    fill = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 13),
    axis.title.x = element_text(size = 13)
  )
