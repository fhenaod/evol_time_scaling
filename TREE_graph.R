library(tidyverse); library(ggthemes)
library(ggpubr); library(cowplot)
library(gg.gap)
library(wesanderson)

# load data ####
path <- c("data2/")
ho <- read.delim(paste0(path,"Ho_et_al_2007_syst_biol.csv"), sep = ",", header = T)
uyeda <- read.delim(paste0(path,"Dryad7.csv"), sep = ",", header = T)
uyeda_2 <- read.delim(paste0(path, "Phylogeniesbynode.csv"), sep = ",", header = T)
henao_diaz <- read.delim(paste0(path, "summary_tree_results.txt"), sep = ",", header = T)
pnas_fossil <- read.delim(paste0(path, "summary_paleo_results.txt"), sep = ",", header = T)

uyeda_3 <- data.frame(rbind(cbind(uyeda$d, uyeda$log10.year), cbind(uyeda_2$d, uyeda_2$Log10Year)))
colnames(uyeda_3) <-  c("d", "log10year")

# individual plots #### 
# molecular
ho_g <- ggplot(ho, aes(x = (cal_kyr/1000), y = (rate_subs.site.myr))) + 
  geom_point(colour = wes_palette("Darjeeling1", n = 3)[1]) + theme_classic() + 
	labs(title = "Molecular", 
	     x = "Calibration time (Myr)", 
	     y = expression(bold(paste("Substution rate (site ", Myr^-1,")"), sep = " "))) + 
	geom_smooth(method = 'nls', formula = y~a*x^b, 
	            method.args = list(start = c(a = .1,b = 1)), se = F, na.rm = T, 
	            #colour="#005fd1", 
	            #colour = wes_palette("Darjeeling1", n = 3)[1],
	            colour = "black",
	            size = 1) + 
	theme_tufte(base_family = "Helvetica") + 
  theme(axis.title.x = element_text(size = 10, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 10, face = "bold", colour = "black"),
        title = element_text(colour = wes_palette("Darjeeling1", n = 3)[1], 
                             face = "bold", size = 12, hjust = .5)) +
  #geom_rangeframe(data = data.frame(x = c(0, 1000), y = c(0, 1)), aes(x, y)) + 
  geom_rangeframe(data = data.frame(x = c(0, 1), y = c(0, 1)), aes(x, y)) + 
  theme(plot.title = element_text(hjust = 0.5))

ho_ins <- ggplot(ho, aes(x = log(cal_kyr/1000), y = log(rate_subs.site.myr))) + 
  geom_point(colour = wes_palette("Darjeeling1", n = 3)[1]) + theme_classic() + 
	labs(title = "", 
	     x = "Log Calibration time (Myr)", y = expression(bold(paste("Log Substution rate (site ", Myr^-1,")"), sep = " "))) + 
  geom_smooth(method = lm, se = F, na.rm = T, 
              #colour = wes_palette("Darjeeling1", n = 3)[1], 
              colour = "black",
              size = 1) + 
	theme_tufte(base_family = "Helvetica") + 
  #geom_rangeframe(data= data.frame(x = c(2, 7), y = c(-6, 0)), aes(x, y)) + 
  geom_rangeframe(data= data.frame(x = c(-5, 0), y = c(-6, 0)), aes(x, y)) + 
  theme(axis.title.x = element_text(size = 10, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 10, face = "bold", colour = "black")) #+
  scale_x_continuous(labels = NULL) + scale_y_continuous(labels = NULL)

# morphology
uyeda_g <- 
  filter(uyeda, Dimensionality..k.==1) %>% 
  filter(Darwins!=0) %>% 
  ggplot(aes(x = (Years/1e6), y = (abs(Darwins)))) + 
            geom_point(colour = wes_palette("Darjeeling1", n = 3)[2]) + theme_classic() + 
									labs(title = "Morphology", 
									     x = "Million years", y = "Absolute rate of evolution (Darwin)") + 
									geom_smooth(method = 'nls', formula = y~a*x^b, 
									            method.args = list(start = c(a = .1, b = -2)), 
									            se = F, na.rm = T, 
									            #colour="#06B95A", 
									            #colour = wes_palette("Darjeeling1", n = 3)[2],
									            colour = "black",
									            size = 1) + 
									theme_tufte(base_family = "Helvetica") + 
  theme(axis.title.x = element_text(size = 10, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 10, face = "bold", colour = "black"),
        title = element_text(colour = wes_palette("Darjeeling1", n = 3)[2], 
                             face = "bold", size = 12, hjust = .5)) +
  #scale_x_continuous(limits = c(0, 1e+05)) + 
  scale_x_continuous(limits = c(0.00001, 100)) + 
  #scale_y_continuous(limits = c(0, 65000)) + 
  #geom_rangeframe(data = data.frame(x = c(0,  1e+05/1000), y = c(0, 65000)), aes(x, y)) + 
  theme(plot.title = element_text(hjust = 0.5))

uyeda_g <- 
  uyeda_g %>% gg.gap(ylim = c(0, 60000), 
                 segments = list(c(1, 1000)), margin = c(0,0,0,.5) 
                 #,tick_width = .1 
                 )
#uyeda_g

uyeda_ins <- filter(uyeda, Dimensionality..k.==1) %>% 
  #mutate(is_more = Darwins >= 1000) %>% 
  filter(Darwins!=0, #is_more == F
         ) %>% 
  ggplot(aes(x = log(Years/1e6), y = log((Darwins)))) + 
  geom_point(colour = wes_palette("Darjeeling1", n = 3)[2]) + 
  theme_classic() + 
  labs(title = "", x = "Log Million years", 
       y = "Log Absolute rate of evolution (Darwin)") + 
  geom_smooth(method = lm, se = F, na.rm = T, 
              #colour="#06B95A", 
              #colour = wes_palette("Darjeeling1", n = 3)[2],
              colour = "black",
              size = 1) + 
									theme_tufte(base_family = "Helvetica") + 
  #scale_x_continuous(limits = c(0, 12.5)) + 
  scale_x_continuous(limits = c(-15, 5)) + 
  scale_y_continuous(limits = c(-8, 15)) +
	geom_rangeframe(data = data.frame(x = c(-15, 5), 
	                                  y = c(-8, 15)), aes(x, y), na.rm = T) +
  theme(axis.title.x = element_text(size = 10, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 10, face = "bold", colour = "black")) #+
  scale_x_continuous(labels = NULL) + scale_y_continuous(labels = NULL)

# phylogenetic
henao_diaz_g <- ggplot(henao_diaz, aes(x = (tree.max.age), y = (mean.clade.lambda))) + 
  geom_point(colour = wes_palette("Darjeeling1", n = 3)[3]) + theme_classic() + 
	labs(title = "Phylogenetic", 
	     x = "Clade age (Myr)", y = expression(bold(paste("Speciation rate (species ", Myr^-1,")"), sep = " "))) + 
	geom_smooth(method = 'nls', formula = y~a*x^b, method.args = list(start = c(a=1,b=2)), se=F, na.rm = T, 
	            #colour="#EA3770", 
	            #colour = wes_palette("Darjeeling1", n = 3)[3],
	            colour = "black",
	            size = 1) + 
	theme_tufte(base_family = "Helvetica") + 
  theme(axis.title.x = element_text(size = 10, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 10, face = "bold", colour = "black"),
        title = element_text(colour = wes_palette("Darjeeling1", n = 3)[3], 
                             face = "bold", size = 12, hjust = .5)) +
  geom_rangeframe(data=data.frame(x = c(0, 350), y = c(0, 1.5)), aes(x, y)) + 
  theme(plot.title = element_text(hjust = 0.5))
	
henao_diaz_ins <- ggplot(henao_diaz, aes(x = log(tree.max.age), y = log(mean.clade.lambda))) + 
  geom_point(colour = wes_palette("Darjeeling1", n = 3)[3]) + theme_classic() + 
	labs(title = "", 
	     x = "Log Clade age (Myr)", y = expression(bold(paste("Log Speciation rate (species ", Myr^-1,")"), sep = " "))) + 
  geom_smooth(method = lm, se = F, na.rm = T, 
              #colour="#EA3770",
              #colour = wes_palette("Darjeeling1", n = 3)[3],
              colour = "black",
              size = 1) + 
	theme_tufte(base_family = "Helvetica") + 
  theme(axis.title.x = element_text(size = 10, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 10, face = "bold", colour = "black")) +
  geom_rangeframe(data=data.frame(x = c(1, 6), y = c(-4, 1)), aes(x, y)) #+ 
  scale_x_continuous(labels = NULL) + scale_y_continuous(labels = NULL)

# fossil 
pnas_g <- ggplot(pnas_fossil, aes(x = (Duration), y = (mean.clade.origination))) + 
  geom_point(colour = wes_palette("Darjeeling2", n = 3)[2]) + 
  theme_classic() + 
  labs(title = "Fossil", 
       x = "Clade duration (Myr)", y = expression(bold(paste("Origination rate (genera ", Myr^-1,")"), sep = " "))) + 
  geom_smooth(method = 'nls', formula = y~a*x^b, 
              method.args = list(start = c(a = .1,b = 1)), se = F, na.rm = T, 
              #colour="#EA3770", 
              #colour = wes_palette("Darjeeling2", n = 3)[2],
              colour = "black",
              size = 1) + 
  theme_tufte(base_family = "Helvetica") + 
  theme(axis.title.x = element_text(size = 10, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 10, face = "bold", colour = "black"),
        title = element_text(colour = wes_palette("Darjeeling2", n = 3)[2], 
                             face = "bold", size = 12, hjust = .5)) +
  geom_rangeframe(data = data.frame(x = c(0, 500), y = c(0, .4)), aes(x, y)) + 
  theme(plot.title = element_text(hjust = 0.5))

pnas_ins <- ggplot(pnas_fossil, aes(x = log(Duration), y = log(mean.clade.origination))) + 
  geom_point(colour = wes_palette("Darjeeling2", n = 3)[2], na.rm = T) + 
  theme_classic() + 
  labs(title = "", 
       x = "Log Clade duration (Myr)", y = expression(bold(paste("Log Origination rate (genera ", Myr^-1,")"), sep = " "))) + 
  geom_smooth(method = lm, se = F, na.rm = T, 
              #colour ="#EA3770", 
              #colour = wes_palette("Darjeeling2", n = 3)[2],
              colour = "black",
              size = 1) + 
  theme_tufte(base_family = "Helvetica") + 
  theme(axis.title.x = element_text(size = 10, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 10, face = "bold", colour = "black")) +
  geom_rangeframe(data = data.frame(x = c(2.5, 6), y = c(-4, -1)), aes(x, y)) + 
  scale_x_continuous(limits = c(2.5, 6)
                     #, labels = NULL
                     ) #+ scale_y_continuous(labels = NULL)

## 8 panel fig ####
gg1 <- ggdraw() + draw_plot(ho_g, x = 0, y = 0) #+ draw_plot(ho_ins, x = .55, y = .6, width = .45, height = .35)
gg2 <- ggdraw() + draw_plot(uyeda_g, x = 0, y = 0) #+ draw_plot(uyeda_ins, x = .55, y = .6, width = .45, height = .35)
gg3 <- ggdraw() + draw_plot(henao_diaz_g, x = 0, y = 0) #+ draw_plot(henao_diaz_ins, x = .55, y = .6, width = .45, height = .35)
gg4 <- ggdraw() + draw_plot(pnas_g, x = 0, y = 0) #+ draw_plot(pna_ins, x = .55, y = .6, width = .45, height = .35)
gg5 <- ggdraw() + draw_plot(ho_ins, x = 0, y = 0) 
gg6 <- ggdraw() + draw_plot(uyeda_ins, x = 0, y = 0)
gg7 <- ggdraw() + draw_plot(henao_diaz_ins, x = 0, y = 0)
gg8 <- ggdraw() + draw_plot(pnas_ins, x = 0, y = 0)

ggarrange(gg1, gg2, gg3, gg4, gg5, gg6, gg7, gg8,
		 labels = LETTERS[1:8],
		 ncol = 4, nrow = 2)
ggsave("8panel_tree.png", dpi = 300, bg = "white",
       width = 35, height = 20, units = "cm",
       limitsize = F)

## 4 panel fig #####
ho_ins <- ggplot(ho, aes(x = log(cal_kyr), y = log(rate_subs.site.myr))) + 
  geom_point(#colour = wes_palette("Darjeeling1", n = 3)[1]
             ) + theme_classic() + 
  labs(title = "", 
       x = "", y = "") + 
  geom_smooth(method = lm, se = F, na.rm = T, 
              colour = wes_palette("Darjeeling1", n = 3)[1], 
              #colour = "black",
              size = 1) + 
  theme_tufte(base_family = "Helvetica") + 
  geom_rangeframe(data= data.frame(x = c(2, 7), y = c(-6, 0)), aes(x, y)) + 
  #theme(axis.title.x = element_text(size = 10, face = "bold", colour = "black"),
  #      axis.title.y = element_text(size = 10, face = "bold", colour = "black")) +
  scale_x_continuous(labels = NULL) + scale_y_continuous(labels = NULL)

uyeda_ins <-  filter(uyeda, Dimensionality..k.==1) %>% 
  filter(Darwins!=0) %>% 
  ggplot(aes(x = log(Years), y = log((Darwins)))) + 
  geom_point(#colour = wes_palette("Darjeeling1", n = 3)[2]
    ) + 
  theme_classic() + 
  labs(title = "", x = "", y = "") + 
  geom_smooth(method = lm, se = F, na.rm = T, 
              #colour="#06B95A", 
              colour = wes_palette("Darjeeling1", n = 3)[2],
              #colour = "black",
              size = 1) + 
  theme_tufte(base_family = "Helvetica") + 
  scale_x_continuous(limits = c(0, 15 )) + scale_y_continuous(limits = c(-8, 15)) +
  geom_rangeframe(data = data.frame(x = c(0, 15), y = c(-8, 15)), aes(x, y), na.rm = T) +
  #theme(axis.title.x = element_text(size = 10, face = "bold", colour = "black"),
  #      axis.title.y = element_text(size = 10, face = "bold", colour = "black")) +
  scale_x_continuous(labels = NULL) + scale_y_continuous(labels = NULL)

henao_diaz_ins <- ggplot(henao_diaz, aes(x = log(tree.max.age), y = log(mean.clade.lambda))) + 
  geom_point(#colour = wes_palette("Darjeeling1", n = 3)[3]
    ) + theme_classic() + 
  labs(title = "", 
       x = "", y = "") + 
  geom_smooth(method = lm, se = F, na.rm = T, 
              #colour="#EA3770",
              colour = wes_palette("Darjeeling1", n = 3)[3],
              #colour = "black",
              size = 1) + 
  theme_tufte(base_family = "Helvetica") + 
  #theme(axis.title.x = element_text(size = 10, face = "bold", colour = "black"),
  #      axis.title.y = element_text(size = 10, face = "bold", colour = "black")) +
  geom_rangeframe(data=data.frame(x = c(1, 6), y = c(-4, 1)), aes(x, y)) + 
  scale_x_continuous(labels = NULL) + scale_y_continuous(labels = NULL)

pnas_ins <- ggplot(pnas_fossil, aes(x = log(Duration), y = log(mean.clade.origination))) + 
  geom_point(#colour = wes_palette("Darjeeling2", n = 3)[2]
    ) + 
  theme_classic() + 
  labs(title = "", 
       x = "", y = "") + 
  geom_smooth(method = lm, se = F, na.rm = T, 
              #colour ="#EA3770", 
              colour = wes_palette("Darjeeling2", n = 3)[2],
              #colour = "black",
              size = 1) + 
  theme_tufte(base_family = "Helvetica") + 
  #theme(axis.title.x = element_text(size = 10, face = "bold", colour = "black"),
  #      axis.title.y = element_text(size = 10, face = "bold", colour = "black")) +
  geom_rangeframe(data = data.frame(x = c(2.5, 5.5), y = c(-4, -1)), aes(x, y), ) + 
  scale_x_continuous(limits = c(2.5, 5.5), labels = NULL) + scale_y_continuous(labels = NULL)

gg5 <- ggdraw() + draw_plot(ho_ins, x = 0, y = 0) 
gg6 <- ggdraw() + draw_plot(uyeda_ins, x = 0, y = 0)
gg7 <- ggdraw() + draw_plot(henao_diaz_ins, x = 0, y = 0)
gg8 <- ggdraw() + draw_plot(pnas_ins, x = 0, y = 0)

gg1 <- ggdraw() + draw_plot(ho_g, x = 0, y = 0) + draw_plot(ho_ins, x = .55, y = .6, width = .45, height = .35)
gg2 <- ggdraw() + draw_plot(uyeda_g, x = 0, y = 0) + draw_plot(uyeda_ins, x = .55, y = .6, width = .45, height = .35)
gg3 <- ggdraw() + draw_plot(henao_diaz_g, x = 0, y = 0) + draw_plot(henao_diaz_ins, x = .55, y = .6, width = .45, height = .35)
gg4 <- ggdraw() + draw_plot(pnas_g, x = 0, y = 0) + draw_plot(pnas_ins, x = .55, y = .6, width = .45, height = .35)

ggarrange(gg1, gg2, gg3, gg4,# gg5, gg6, gg7, gg8,
          labels = LETTERS[1:4],
          ncol = 4, nrow = 1)
ggsave("4panel_tree.png", dpi = 300, bg = "white",
       width = 45, height = 20, units = "cm",
       limitsize = F)


# log-log transformed 

ho_g <- ggplot(ho, aes(x = log(cal_kyr), y = log(rate_subs.site.myr))) + geom_point() + theme_classic() + 
	labs(subtitle = "Molecular", x = "Ln Calibration time (kyr)", y = expression(paste("Ln Substution rate (site ", Myr^-1,")"), sep = " ")) + 
	geom_smooth(method = 'nls', formula = y~a*x^b, method.args = list(start = c(a = 1,b = 1)), se = F, na.rm = T, colour="#005fd1", size = 1.5) + 
	theme_tufte(base_family = "Helvetica") + geom_rangeframe(data=data.frame(x = c(1, 7), y = c(-5, 0)), aes(x, y))

uyeda_g <-  ggplot(uyeda, aes(x = log(Years), y = log(abs(Darwins)))) + geom_point() + theme_classic() + 
									labs(subtitle = "Morphology", x = "Ln Years", y = "Ln Rate of evolution (Darwin)") + 
									geom_smooth(method = 'nls', formula = y~a*x^b, method.args = list(start = c(a = 1,b = 1)), se = F, na.rm = T, colour="#06B95A", size = 1.5) + 
									theme_tufte(base_family = "Helvetica") + scale_x_continuous(limits = c(0, 15 )) + scale_y_continuous( limits = c(-6, 15)) +
									geom_rangeframe(data = data.frame(x = c(0, 15), y = c(-6, 15)), aes(x, y))

henao_diaz_g <- ggplot(henao_diaz, aes(x = log(tree.max.age), y = log(mean.clade.lambda))) + geom_point() + theme_classic() + 
	labs(subtitle = "Phylogenetic", x = "Ln Clade age (Myr)", y = expression(paste("Ln Speciation rate (species ", Myr^-1,")"), sep = " ")) + 
	geom_smooth(method = 'nls', formula = y~a*x^b, method.args = list(start = c(a=1,b=2)), se=F, na.rm = T, colour="#EA3770", size=1.5) + 
	theme_tufte(base_family = "Helvetica") + geom_rangeframe(data=data.frame(x = c(1, 6), y = c(-4, 1)), aes(x, y))

ggarrange(ho_g, uyeda_g, henao_diaz_g,
		 labels = c("A", "B", "C"),
		 ncol = 3, nrow = 1)
		 	 
### uyeda's d plots
		
uyeda_gl1<- ggplot(uyeda_3, aes(x = log10year, y = d)) + geom_point() + theme_classic() + 
									labs(subtitle = "Morphology", x = "log10year", y = "d") + 
									geom_smooth(method = 'nls', formula = y~a*x^b, method.args = list(start = c(a = .5,b = 1)), se = F, na.rm = T, colour="#06B95A", size = 1.5) + 
									theme_tufte(base_family = "Helvetica") 
									#scale_y_continuous( limits = c(0,8)) + scale_x_continuous(limits = c(0,2)) + 
									#geom_rangeframe(data = data.frame(x = c(0, 15), y = c(-50, 50)), aes(x, y))
									
uyeda_gl2 <- ggplot(uyeda_3, aes(x = log10year, y = abs(d/log10year))) + geom_point() + theme_classic() + 
									labs(subtitle = "Morphology", x = "log10year", y = "Ln rate of evolution (Darwin)") + 
									geom_smooth(method = 'nls', formula = y~a*x^b, method.args = list(start = c(a = .5,b = 1)), se = F, na.rm = T, colour="#06B95A", size = 1.5) + 
									theme_tufte(base_family = "Helvetica")  
									#scale_y_continuous( limits = c(0,8)) + scale_x_continuous(limits = c(0,2)) + 
									#geom_rangeframe(data = data.frame(x = c(0, 15), y = c(-50, 50)), aes(x, y))
									
uyeda_gl3 <- ggplot(uyeda_3, aes(x = log10year, y = (d/log10year))) + geom_point() + theme_classic() + 
									labs(subtitle = "Morphology", x = "log10year", y = "Ln rate of evolution (Darwin)") + 
									geom_smooth(method = 'nls', formula = y~a*x^b, method.args = list(start = c(a = .5,b = 1)), se = F, na.rm = T, colour="#06B95A", size = 1.5) + 
									theme_tufte(base_family = "Helvetica") 
									#scale_y_continuous( limits = c(0,8)) + scale_x_continuous(limits = c(0,2)) + 
									#geom_rangeframe(data = data.frame(x = c(0, 15), y = c(-50, 50)), aes(x, y))									
							
ggarrange(uyeda_gl1, uyeda_gl2, uyeda_gl3,
		 labels = c("A", "B", "C"),
		 ncol = 3, nrow = 1)	

# descriptive regressions

lm(log(rate_subs.site.myr)~log(cal_kyr), data = ho) %>% summary()
lm(log(rate_subs.site.myr)~log(cal_kyr), data = ho) %>% confint(level = .95) %>% round(3)

lm(log(mean.clade.lambda)~log(tree.max.age), data = henao_diaz) %>% summary()
lm(log(mean.clade.lambda)~log(tree.max.age), data = henao_diaz) %>% confint(level = .95) %>% round(3)

filter(uyeda, Dimensionality..k.==1) %>% lm(log((Darwins))~log(Years), data = .) %>% summary() %>% coefficients()
