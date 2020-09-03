
require(tidyverse)
require(ggbeeswarm)
require(BuenColors)
require(gtools)
require(lubridate)
library(sjPlot)
library(sjmisc)

#colors
hbb3.5 <- c("#ede1ee", "#d4bad8", "#b688bd")
hbd3.5 <-c("#76d6ff", "#ffc900", "#ff1414")
hbbhbd <- c("#d3f0fc", "#a8e1f9", "#7CD2F6")
bcl11a <- c("#838383", "#ff9300", "#009051", "#ff1414")
            

# read in data
df <- read_csv("~/Desktop/Yongs_modelling/recoded_2.csv") %>%
   mutate(Gamma = G1 + G2) %>%
   mutate(class = paste(BCL11A, HBB_HBD, three.five,Gamma, sep = "_")) %>%
   mutate(Donor = factor(Donor))


# update df with a code
df <- df %>%
      mutate(code = paste0(BCL11A, HBB_HBD, three.five, G1, G2))

# look at donor effects in uneditted
uneditted <- df %>%
   filter(code == "00000")

summary(lm(HBG ~ Donor,
           data = uneditted))

df <- df %>%
   mutate(Donor = ifelse(!(Donor %in% c())))

# only donor 3684, 3734, 3744 have significant effect
# will recode all other donors to be the same (too many donors will result in
# singular fit for lmm)

df <- df %>%
   mutate(Donor = as.character(Donor)) %>%
   mutate(Donor = ifelse(!(Donor %in% c())))

# beta locus
beta_locus_data <- df %>% filter(BCL11A == 0 & Gamma  == 0) %>%
   filter(HBG13bp_gRNA ==0 & BCL11A_exon2_gRNA == 0 & BCL11A_exon4_gRNA == 0)

lm_HBG <- lm(HBG ~ HBB_HBD + three.five + HBB_HBD:three.five, data = beta_locus_data)
lm_HBG <- lmerTest::lmer(HBG ~ HBB_HBD + three.five + HBB_HBD:three.five +
                            (1 | `HBD-3.5kb_gRNA`) , data = beta_locus_data)
summary(lm_HBG)



p <- plot_model(lm_HBG, type = "pred", terms = c("HBB_HBD", "three.five"), colors = hbd3.5) + 
      ylim(0, 1.5) +
      scale_x_continuous(breaks=seq(0,2,1)) +
      pretty_plot()
ggsave("~/Desktop/Yongs_modelling/HBG_HBBHBD_3.5kb.pdf", p, height = 4, width = 5, units = "in")


lm_HBB <- lm(HBB ~ HBB_HBD + three.five + HBB_HBD:three.five, data = temp)
lm_HBB <- lmerTest::lmer(HBB ~ HBB_HBD + three.five + HBB_HBD:three.five +
                  (1 | HBB_gRNA) + (1 |HBD_gRNA) + (1 | `HBD-3.5kb_gRNA`) + (1 | Donor), data = temp)
summary(lm_HBB)
p <- plot_model(lm_HBB, type = "pred", terms = c("HBB_HBD", "three.five"), colors = hbd3.5) + 
      scale_x_continuous(breaks=seq(0,2,1)) +
      pretty_plot()
ggsave("~/Desktop/Yongs_modelling/HBB_HBBHBD_3.5kb.pdf", p, height = 4, width = 5, units = "in")



df2 <- df %>%
      filter(HBB_HBD == three.five) %>%
      mutate(HBB_3.5 = HBB_HBD)

df2 <- df2 %>%
      filter(HBB_HBD == three.five) %>%
      mutate(HBB_3.5 = HBB_HBD) %>%
      mutate(BCL11A =  case_when(grepl("exon4", ID) & BCL11A == 1 ~ 2,
                                 grepl("exon4", ID) & BCL11A == 2 ~ 3,
                                 TRUE ~ BCL11A)) %>%
   filter(HBD_g ==0)

lm_HBG_BCL11A <- lm(HBG ~BCL11A * (HBB_3.5) , data = df2)
lm_HBG <- lmerTest::lmer(HBG ~ BCL11A * HBB_3.5 + (1 | AAVS1_gRNA) + (1 | BCL11A_exon2_gRNA) + (1 | BCL11A_exon4_gRNA)  +
                  (1 | HBB_gRNA) + (1 | `HBD-3.5kb_gRNA`) + (1 | HBG13bp_gRNA) + (1 | Donor), data = df2)
p <- plot_model(lm_HBG, type = "pred", terms = c( "HBB_3.5", "BCL11A"),colors = bcl11a) +
      ylim(0, 1.75) +
      scale_x_continuous(breaks=seq(0,2,1)) +
      pretty_plot()

ggsave("~/Desktop/Yongs_modelling/HBG_HBB35kb_BCl11A.pdf", p, height = 4, width = 5, units = "in")

p <- plot_model(lm_HBG, type = "pred", terms = c( "Gamma", "BCL11A"),colors = bcl11a) +
   ylim(0, 1.75) +
   scale_x_continuous(breaks=seq(0,2,1)) +
   pretty_plot()

ggsave("~/Desktop/Yongs_modelling/HBHBG13bp_gRNAamma_BCl11A.pdf", p, height = 4, width = 5, units = "in")


lm_HBB <- lmer(HBB  ~ BCL11A * (HBB_3.5 + Gamma) + (1 | AAVS1_gRNA) + (1 | BCL11A_exon2_gRNA) + (1 | BCL11A_exon4_gRNA)  +
                  (1 | HBB_gRNA) + (1 | `HBD-3.5kb_gRNA`) + (1 | HBG13bp_gRNA) + (1 | Donor), data = df2)
p <- plot_model(lm_HBB, type = "pred", terms = c( "HBB_3.5", "BCL11A"),colors = bcl11a) +
      scale_x_continuous(breaks=seq(0,2,1)) +
      pretty_plot()

ggsave("~/Desktop/Yongs_modelling/HBB_HBB35kb_BCl11A.pdf", p, height = 4, width = 5, units = "in")


p <- plot_model(lm_HBB, type = "pred", terms = c( "Gamma", "BCL11A"),colors = bcl11a) +
   scale_x_continuous(breaks=seq(0,4,1)) +
   pretty_plot()

ggsave("~/Desktop/Yongs_modelling/HBB_gRNAamma_BCl11A.pdf", p, height = 4, width = 5, units = "in")

lm_HBB <- lmer(HBB ~ BCL11A * (HBB_3.5 + Gamma) + (1 | AAVS1_gRNA) + (1 | BCL11A_exon2_gRNA) + (1 | BCL11A_exon4_gRNA)  +
                  (1 | HBB_gRNA) + (1 | `HBD-3.5kb_gRNA`) + (1 | HBG13bp_gRNA) + (1 | Batch), data = df2)
p <- plot_model(lm_HBB, type = "pred", terms = c( "Gamma", "BCL11A"),colors = bcl11a) +
      scale_x_continuous(breaks=seq(0,4,1)) +
      pretty_plot()

ggsave("~/Desktop/Yongs_modelling/HBB_gRNAamma_BCl11A.pdf", p, height = 4, width = 5, units = "in")


lm_HBB_BCL11A <- lm(HBB ~BCL11A * (Gamma) , data = df2)
p <- plot_model(lm_HBB_BCL11A, type = "pred", terms = c( "Gamma", "BCL11A"),colors = bcl11a) +
      ylim(-0.25, 1) +
      scale_x_continuous(breaks=seq(0,4,1)) +
      pretty_plot()

ggsave("~/Desktop/Yongs_modelling/HBB_gRNAamma_BCl11A.pdf", p, height = 4, width = 5, units = "in")

lm(HBG ~BCL11A * (HBB_3.5 + Gamma) , data = df2)

int_hbb <- expand.grid(0:2, 0:3) %>%
      rename_all(~c("HBB_3.5", "BCL11A")) %>%
      mutate(Intercept = 0.1672286,
             n_interactions = HBB_3.5 * BCL11A,
             HBB_cont = HBB_3.5 * 0.4210317,
             BCL11A_cont = BCL11A *  0.4039013,
             Int_cont = -0.1450618*HBB_3.5*BCL11A,
             predicted = HBB_cont + BCL11A_cont + Int_cont)


int_gamma <- expand.grid(0:4, 0:3) %>%
      rename_all(~c("Gamma", "BCL11A")) %>%
      mutate(Intercept = 0.1672286,
             n_interactions = Gamma * BCL11A,
             Gamma_cont = Gamma* 0.2064,
             BCL11A_cont = BCL11A *  0.4039013,
             Int_cont = -0.1392*Gamma*BCL11A,
             predicted = Gamma_cont + BCL11A_cont + Int_cont)



