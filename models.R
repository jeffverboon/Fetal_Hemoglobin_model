
require(tidyverse)
require(ggbeeswarm)
require(BuenColors)
require(gtools)
require(lubridate)
library(sjPlot)
library(sjmisc)
library(lme4)

#colors
hbb3.5 <- c("#ede1ee", "#d4bad8", "#b688bd")
hbd3.5 <-c("#76d6ff", "#ffc900", "#ff1414")
hbbhbd <- c("#d3f0fc", "#a8e1f9", "#7CD2F6")
bcl11a <- c("#BBBBBB", "#F2B669", "#59BA8C", "#E82B63")
zbtb7a <- c("#AAAAAA", "#008DFF", "#EC67FF")


# read in data
df <- read_csv("./data.csv") %>%
   mutate(Donor = factor(Donor))


# update df with a code for edits
df <- df %>%
   mutate(code = paste0(ZBTB7A,BCL11A, HBB_HBD, three.five, G1, G2))

# look at donor effects in uneditted
uneditted <- df %>%
   filter(code == "000000")

summary(lm(HBG ~ Donor,
           data = uneditted))
summary(lm(HBB ~ Donor,
           data = uneditted))



# only donor 3684, 3734, 3744 have significant effect
# will recode all other donors to be the same (too many donors will result in
# singular fit for lmm)

df <- df %>%
   mutate(Donor = as.character(Donor)) %>%
   mutate(Donor = factor(ifelse(!(Donor %in% c("3684", "3734","3744")), "Normal", Donor)))

# lets see if the gRNA's are affecting the baseline hbf. 
# we can test this by seeing if in uneditted conditions
# changes in the HBB or HBG. 

summary(lm(HBG ~ AAVS1_gRNA + BCL11A_exon2_gRNA + BCL11A_exon4_gRNA + HBB_gRNA +
              HBD_gRNA + `HBD-3.5kb_gRNA` + HBG13bp_gRNA + ZBTB7A_gRNA,
           data = uneditted))
summary(lm(HBB ~ AAVS1_gRNA + BCL11A_exon2_gRNA + BCL11A_exon4_gRNA + HBB_gRNA +
              HBD_gRNA + `HBD-3.5kb_gRNA` + HBG13bp_gRNA + ZBTB7A_gRNA,
           data = uneditted))

# there does appear to be some affect with signicance dependant on numbers
# good candidate for mixed model of the intercept


# beta locus
beta_locus_data <- df %>% filter(BCL11A == 0 & G1  == 0 & G2 == 0 & ZBTB7A ==0) %>%
   filter(HBG13bp_gRNA ==0 & BCL11A_exon2_gRNA == 0 & BCL11A_exon4_gRNA == 0)

# fit model
lm_HBG <- lmerTest::lmer(HBG ~ HBB_HBD + three.five + HBB_HBD:three.five +
                            (1 | HBB_gRNA) + (1 |HBD_gRNA) + (1 | AAVS1_gRNA) +
                            (1 | `HBD-3.5kb_gRNA`) + (1 | Donor), data = beta_locus_data)
summary(lm_HBG)


# determine whcih random effects are causing singular fit and remove from model
summary(rePCA(lm_HBG))
lm_HBG <- lmerTest::lmer(HBG ~ HBB_HBD + three.five + HBB_HBD:three.five +
                            (1 | `HBD-3.5kb_gRNA`) + (1 | Donor), data = beta_locus_data)

p <- plot_model(lm_HBG, type = "pred", terms = c("HBB_HBD", "three.five"), colors = hbd3.5) + 
   ylim(0, 1.5) +
   scale_x_continuous(breaks=seq(0,2,1)) +
   pretty_plot()
ggsave("./HBG_HBBHBD_3.5kb.pdf", p, height = 4, width = 5, units = "in")


# the HBB-HBD condition is removing copies of adult globin; model as well
lm_HBB <- lmerTest::lmer(HBB ~ HBB_HBD + three.five + HBB_HBD:three.five +
                            (1 | HBB_gRNA) + (1 |HBD_gRNA) + (1 | AAVS1_gRNA) +
                            (1 | `HBD-3.5kb_gRNA`) + (1 | Donor), data = beta_locus_data)
# fix singular fit
summary(rePCA(lm_HBB))
lm_HBB <- lmerTest::lmer(HBB ~ HBB_HBD + three.five + HBB_HBD:three.five +
                            (1 | HBB_gRNA) + (1 | `HBD-3.5kb_gRNA`) +
                            (1 | Donor), data = beta_locus_data)
summary(lm_HBB)
p <- plot_model(lm_HBB, type = "pred", terms = c("HBB_HBD", "three.five"), colors = hbd3.5) + 
   scale_x_continuous(breaks=seq(0,2,1)) +
   pretty_plot()
ggsave("./HBB_HBBHBD_3.5kb.pdf", p, height = 4, width = 5, units = "in")


# for edits combined with BCL11A or Gamma we only did the HBB-3.5kb deletion; 
# remove partial deletions from dataframe and make new variable for that
# remove uneditted that were were treated with HBD_gRNA as well
# combined model
# remove new ZBTB7A a well for this model

df2 <- df %>%
   filter(HBB_HBD == three.five) %>%
   mutate(HBB_3.5 = HBB_HBD) %>%
   filter(HBD_gRNA ==0 & ZBTB7A == 0) 

# we also will treat gamma1 and gamma2 as a combined doseage (0:4)
# and BCL11A as 0:3 as described in the manuscript

df2 <- df2 %>%
   mutate(Gamma = G1 + G2) %>%
   mutate(BCL11A =  case_when(grepl("exon4", ID) & BCL11A == 1 ~ 2,
                              grepl("exon4", ID) & BCL11A == 2 ~ 3,
                              TRUE ~ BCL11A)) 

#
lm_HBG <- lmerTest::lmer(HBG ~ BCL11A * (Gamma) + 
                  (1 | HBG13bp_gRNA), data = df2a)


# model HbF and BCL11A + adult beta region altogther
lm_HBG <- lmerTest::lmer(HBG ~ BCL11A * (Gamma + HBB_3.5) + (1 | AAVS1_gRNA) + (1 | BCL11A_exon2_gRNA) + (1 | BCL11A_exon4_gRNA)  +
                            (1 | HBB_gRNA) + (1 | `HBD-3.5kb_gRNA`) + (1 | HBG13bp_gRNA) + (1 | Donor), data = df2)

# fix singular fit
summary(rePCA(lm_HBG))
lm_HBG <- lmerTest::lmer(HBG ~ BCL11A * (Gamma + HBB_3.5) + (1 | BCL11A_exon4_gRNA) + 
                            (1 | HBG13bp_gRNA), data = df2)


p <- plot_model(lm_HBG, type = "pred", terms = c( "HBB_3.5", "BCL11A"),colors = bcl11a) +
   ylim(0, 1.75) +
   scale_x_continuous(breaks=seq(0,2,1)) +
   pretty_plot()

ggsave("./HBG_HBB35kb_BCl11A.pdf", p, height = 4, width = 5, units = "in")

p <- plot_model(lm_HBG, type = "pred", terms = c( "Gamma", "BCL11A"),colors = bcl11a) +
   ylim(0, 1.75) +
   scale_x_continuous(breaks=seq(0,4,1)) +
   pretty_plot()

ggsave("./HBG_gamma_BCl11A.pdf", p, height = 4, width = 5, units = "in")


lm_HBB <- lmerTest::lmer(HBB  ~ BCL11A * (Gamma + HBB_3.5) + (1 | AAVS1_gRNA) + (1 | BCL11A_exon2_gRNA) + (1 | BCL11A_exon4_gRNA)  +
                            (1 | HBB_gRNA) + (1 | `HBD-3.5kb_gRNA`) + (1 | HBG13bp_gRNA) + (1 | Donor), data = df2)

# singularity fix
summary(rePCA(lm_HBB))
lm_HBB <- lmerTest::lmer(HBB  ~ BCL11A * (Gamma + HBB_3.5) + (1 | AAVS1_gRNA) + (1 | BCL11A_exon4_gRNA)  +
                             (1 | Donor), data = df2)


p <- plot_model(lm_HBB, type = "pred", terms = c( "HBB_3.5", "BCL11A"),colors = bcl11a) +
   scale_x_continuous(breaks=seq(0,2,1)) +
   pretty_plot()

ggsave("./HBB_HBB35kb_BCl11A.pdf", p, height = 4, width = 5, units = "in")


p <- plot_model(lm_HBB, type = "pred", terms = c( "Gamma", "BCL11A"),colors = bcl11a) +
   scale_x_continuous(breaks=seq(0,4,1)) +
   pretty_plot()

ggsave("./HBB_gamma_BCl11A.pdf", p, height = 4, width = 5, units = "in")

# ZBTB7A now

df3 <- df %>%
   mutate(Gamma = G1 + G2) %>%
   filter(HBB_HBD == three.five) %>%
   mutate(HBB_3.5 = HBB_HBD) %>%
   filter(Gamma == 0) %>%
   mutate(BCL11A =  case_when(grepl("exon4", ID) & BCL11A == 1 ~ 2,
                              grepl("exon4", ID) & BCL11A == 2 ~ 3,
                              TRUE ~ BCL11A))

# first model HBG
# need to add in the interaction between BCL11A and HBB_3.5 as well
lm_HBG <- lmer(HBG ~ ZBTB7A + BCL11A + HBB_3.5 + ZBTB7A:BCL11A +
                  ZBTB7A:HBB_3.5 + BCL11A:HBB_3.5 + (1 | AAVS1_gRNA) + 
                  (1 | BCL11A_exon2_gRNA) + (1 | BCL11A_exon4_gRNA)  +
                  (1 | HBB_gRNA) + (1 | `HBD-3.5kb_gRNA`) + (1 | HBG13bp_gRNA) +
                  (1 | Donor) + (1 | ZBTB7A_gRNA), data = df3)

# fix singular fit
summary(rePCA(lm_HBG))
lm_HBG <- lmerTest::lmer(HBG ~ ZBTB7A + BCL11A + HBB_3.5 + ZBTB7A:BCL11A +
                            ZBTB7A:HBB_3.5 + BCL11A:HBB_3.5 + (1 | BCL11A_exon4_gRNA),
                         data = df3)


p <- plot_model(lm_HBG, type = "pred", terms = c("BCL11A", "ZBTB7A"), colors = zbtb7a) +
   scale_x_continuous(breaks=seq(0,2,1)) +
   xlim(0,2) +
   pretty_plot()

ggsave("./HBG_ZBTB7A_BCl11A.pdf", p, height = 4, width = 5, units = "in")

p <- plot_model(lm_HBG, type = "pred", terms = c("HBB_3.5", "ZBTB7A"), colors = zbtb7a) +
   scale_x_continuous(breaks=seq(0,2,1)) +
   pretty_plot()

ggsave("./HBG_ZBTB7A_HBB35kb.pdf", p, height = 4, width = 5, units = "in")


lm_HBB <- lmerTest::lmer(HBB  ~ ZBTB7A + BCL11A + HBB_3.5 + ZBTB7A:BCL11A +
                            ZBTB7A:HBB_3.5 + BCL11A:HBB_3.5 + (1 | AAVS1_gRNA) + 
                            (1 | BCL11A_exon2_gRNA) + (1 | BCL11A_exon4_gRNA)  +
                            (1 | HBB_gRNA) + (1 | `HBD-3.5kb_gRNA`) + (1 | HBG13bp_gRNA) +
                            (1 | Donor) + (1 | ZBTB7A_gRNA), data = df3)
# singularity fix
summary(rePCA(lm_HBB))
lm_HBB <- lmerTest::lmer(HBB  ~ ZBTB7A + BCL11A + HBB_3.5 + ZBTB7A:BCL11A +
                            ZBTB7A:HBB_3.5 + BCL11A:HBB_3.5 + (1 | AAVS1_gRNA) + 
                            (1 | BCL11A_exon4_gRNA)  +
                            (1 | HBB_gRNA) +
                            (1 | Donor), data = df3)

p <- plot_model(lm_HBB, type = "pred", terms = c("BCL11A","ZBTB7A"), colors = zbtb7a) +
   scale_x_continuous(breaks=seq(0,2,1)) +
   xlim(0,2) +
   pretty_plot()

ggsave("./HBB_ZBTB7A_BCl11A.pdf", p, height = 4, width = 5, units = "in")


p <- plot_model(lm_HBB, type = "pred", terms = c("HBB_3.5", "ZBTB7A"), colors = zbtb7a) +
   scale_x_continuous(breaks=seq(0,2,1)) +
   pretty_plot()

ggsave("./HBB_ZBTB7A_HBB3.5.pdf", p, height = 4, width = 5, units = "in")
