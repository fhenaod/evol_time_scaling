library(tidyverse)
library(ggthemes)
library(ggpubr)
library(cowplot)

path <- c("~/Dropbox/Proyectos/UBC/Biology/Harmon_et_al_AREEs/data/")
ho <- read.delim(paste0(path,"Ho_et_al_2007_syst_biol.csv"), sep = ",", header = T)
uyeda <- read.delim(paste0(path,"Dryad7.csv"), sep = ",", header = T)
uyeda_2 <- read.delim(paste0(path, "Phylogeniesbynode.csv"), sep = ",", header = T)
henao_diaz <- read.delim(paste0(path, "summary_tree_results.txt"), sep = ",", header = T)

uyeda_3 <- data.frame(rbind(cbind(uyeda$d, uyeda$log10.year), cbind(uyeda_2$d, uyeda_2$Log10Year)))
colnames(uyeda_3) <-  c("d", "log10year")

# final AREEs paper 
ho_g <- ggplot(ho, aes(x = (cal_kyr), y = (rate_subs.site.myr))) + geom_point() + theme_classic() + 
	labs(subtitle = "Molecular", x = "Calibration time (kyr)", y = expression(paste("Substution rate (site ", Myr^-1,")"), sep = " ")) + 
	geom_smooth(method = 'nls', formula = y~a*x^b, method.args = list(start = c(a = 1,b = 1)), se = F, na.rm = T, colour="#005fd1", size = 1.5) + 
	theme_tufte(base_family = "Helvetica") + geom_rangeframe(data=data.frame(x = c(0, 1000), y = c(0, 1)), aes(x, y))
	
ho_ins <- ggplot(ho, aes(x = log(cal_kyr), y = log(rate_subs.site.myr))) + geom_point() + theme_classic() + 
	labs(subtitle = "", x = "", y = "") + geom_smooth(method = lm, se = F, na.rm = T, colour="#005fd1", size = 1.5) + 
	theme_tufte(base_family = "Helvetica") + geom_rangeframe(data=data.frame(x = c(2, 7), y = c(-6, 0)), aes(x, y)) + scale_x_continuous(labels = NULL) + scale_y_continuous(labels = NULL)

uyeda_g <- filter(uyeda, Dimensionality..k.==1) %>% ggplot(aes(x = (Years), y = (abs(Darwins)))) + geom_point() + theme_classic() + 
									labs(subtitle = "Morphology", x = "Year", y = "Abs Rate of evolution (Darwin)") + 
									geom_smooth(method = 'nls', formula = y~a*x^b, method.args = list(start = c(a = .1,b = 1)), se = F, na.rm = T, colour="#06B95A", size = 1.5) + 
									theme_tufte(base_family = "Helvetica") + scale_y_continuous( limits = c(0,65000)) + scale_x_continuous(limits = c(0,750)) + 
									geom_rangeframe(data = data.frame(x = c(0, 750), y = c(0, 65000)), aes(x, y))
									
uyeda_ins <-  filter(uyeda, Dimensionality..k.==1) %>% ggplot(aes(x = log(Years), y = log((Darwins)))) + geom_point() + theme_classic() + 
									labs(subtitle = "", x = "", y = "") + geom_smooth(method = lm, se = F, na.rm = T, colour="#06B95A", size = 1.5) + 
									theme_tufte(base_family = "Helvetica") + scale_x_continuous(limits = c(0, 15 )) + scale_y_continuous( limits = c(-8, 15)) +
									geom_rangeframe(data = data.frame(x = c(0, 15), y = c(-8, 15)), aes(x, y), na.rm = T) + scale_x_continuous(labels = NULL) + scale_y_continuous(labels = NULL)

henao_diaz_g <- ggplot(henao_diaz, aes(x =(tree.max.age), y =(mean.clade.lambda))) + geom_point() + theme_classic() + 
	labs(subtitle = "Phylogenetic", x = " Clade age (Myr)", y = expression(paste("Speciation rate (species ", Myr^-1,")"), sep = " ")) + 
	geom_smooth(method = 'nls', formula = y~a*x^b, method.args = list(start = c(a=1,b=2)), se=F, na.rm = T, colour="#EA3770", size=1.5) + 
	theme_tufte(base_family = "Helvetica") + geom_rangeframe(data=data.frame(x = c(0, 350), y = c(0, 1.5)), aes(x, y))
	
henao_diaz_ins <- ggplot(henao_diaz, aes(x = log(tree.max.age), y = log(mean.clade.lambda))) + geom_point() + theme_classic() + 
	labs(subtitle = "", x = " ", y = "") + geom_smooth(method = lm, se=F, na.rm = T, colour="#EA3770", size=1.5) + 
	theme_tufte(base_family = "Helvetica") + geom_rangeframe(data=data.frame(x = c(1, 6), y = c(-4, 1)), aes(x, y)) + scale_x_continuous(labels = NULL) + scale_y_continuous(labels = NULL)

gg1 <- ggdraw() + draw_plot(ho_g, x = 0, y = 0) + draw_plot(ho_ins, x = .55, y = .6, width = .45, height = .35)
gg2 <- ggdraw() + draw_plot(uyeda_g, x = 0, y = 0) + draw_plot(uyeda_ins, x = .55, y = .6, width = .45, height = .35)
gg3 <- ggdraw() + draw_plot(henao_diaz_g, x = 0, y = 0) + draw_plot(henao_diaz_ins, x = .55, y = .6, width = .45, height = .35)

ggarrange(gg1, gg2, gg3,
		 labels = c("A", "B", "C"),
		 ncol = 3, nrow = 1)

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
