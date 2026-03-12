# ============================================================
# AUTO-GENERATED — DO NOT EDIT DIRECTLY
# Edits will be overwritten on next org-babel tangle.
# 
# Source:  /home/jeszyman/repos/frag/frag.org
# Author:  Jeffrey Szymanski
# Tangled: 2026-03-10 18:04:56
# ============================================================

source("~/repos/mpnst-frag/config/library_loads.R")
library(tidyverse)

frag_count = read.table("/mnt/ris/aadel/mpnst/frag/frag_counts.tsv", header = F)
load("/mnt/ris/aadel/mpnst/data_model/data_model.RData")

names(frag_count) = c("library_id","frag_length","chr","start","end","count")

test =
  frag_count %>%
  pivot_wider(names_from = frag_length, values_from = count) %>%
  group_by(library_id,chr,start,end) %>%
  mutate(ratio = short/long)

washout_libs = c("lib218","lib107","lib117","lib126","lib129","lib142","lib158","lib175","lib182","lib184","lib202","lib205")


test2 = libraries_full %>%
  filter(library_type == "wgs") %>%
  filter(isolation_type == "cfdna") %>%
  filter(institution %in% c("nci","washu")) %>%
  filter(current_dx %in% c("plexiform","healthy") | library_id %in% washout_libs)


test2 = libraries_full %>%
  filter(library_type == "wgs") %>%
  filter(isolation_type == "cfdna") %>%
  filter(institution %in% c("nci","washu")) %>%
  filter(current_dx %in% c("healthy", "plexiform"))

dx = test2 %>% select(library_id, current_dx)

frags =
  test %>% filter(library_id %in% dx$library_id)

test = frags %>% select(library_id, chr, start, end, ratio) %>% pivot_wider(names_from = library_id, values_from = ratio)

test2 = test
head(test2)
test2[4:91] = scale(test2[4:91])


test3 = test2 %>% pivot_longer(starts_with("lib"), names_to = "library_id", values_to = "ratio") %>% left_join(dx, by = "library_id")

test3 %>% filter(chr == "chr1") %>% ggplot(., aes(x = start, y = ratio, color = current_dx, group = library_id)) +
  geom_line(stat = "smooth", span = 0.1, alpha = 0.8, aes(size = current_dx)) + facet_grid(~chr) + scale_size_manual(values = c(5,.5,.5))


plot =
test3 %>% mutate(new_id = library_id) %>%
mutate(new_id = ifelse(current_dx == "healthy", "healthy", library_id )) %>%
ggplot(., aes(x = start, y = ratio, group = library_id, color = current_dx, linetype = current_dx)) +
  geom_line(stat = "smooth", alpha = 0.8, span = 0.3) + facet_wrap(~chr, ncol = 2, scales = "free") + scale_size_manual(values = c(1,.5,.5))
ggsave(plot, width = 30, height = 40, filename = "/tmp/plot.pdf")


plot2 =
test3 %>% mutate(new_id = library_id) %>%
mutate(new_id = ifelse(current_dx == "healthy", "healthy", library_id )) %>%
ggplot(., aes(x = start, y = ratio, group = current_dx, color = current_dx, linetype = current_dx)) +
  geom_smooth(alpha = 0.8, span = 0.3, aes(fill = current_dx)) + facet_wrap(~chr, ncol = 2, scales = "free")
ggsave(plot2, width = 30, height = 40, filename = "/tmp/plot2.pdf")



 geom_line(stat="smooth",method = "lm", formula = y ~ 0 + I(1/x) + I((x-1)/x),
              size = 1.5,
              linetype ="dashed",
              alpha = 0.5)

test3 %>% filter(chr %in% c("chr20","chr17")) %>% ggplot(., aes(x = start, y = ratio, color = current_dx)) + geom_smooth(se = F, span = .2, alpha = 0.1) + facet_grid(~chr)


head(test3)

head(test2)

mat = test2[,-c(1,2,3)]

mat = as.matrix(mat)

rownames(mat) = paste(test2$chr,test2$start,test2$end,sep = "_")
head(mat)

mat = t(mat)

pca = prcomp(mat)

# Get principle component 1 & 2 values
(pve_pc1=round(100*summary(pca)$importance[2,1]))
(pve_pc2=round(100*summary(pca)$importance[2,2]))

summary(pca)$importance

head(pca$x)

pca_plot = as.data.frame(pca$x) %>%
  rownames_to_column(var = "library_id") %>%
  left_join(dx, by = "library_id") %>%
  ggplot(., aes(x = PC1, y = PC2, color = current_dx)) +
  geom_point(size = 4)
pca_plot

+
  theme_cowplot() +
  xlab(paste("PC1, ", pve_pc1, "% variance explained", sep ="")) +
  ylab(paste("PC2, ", pve_pc2, "% variance explained", sep =""))
pca_plot


pca_plot = as.data.frame(pca$x) %>%
  rownames_to_column(var = "sample_id") %>%
  mutate(cohort_id = ifelse(grepl("a", sample_id), "ir", "sham")) %>%
  ggplot(., aes(x = PC1, y = PC2, color = cohort_id)) +
  geom_point(size = 4) +
  theme_cowplot() +
  xlab(paste("PC1, ", pve_pc1, "% variance explained", sep ="")) +
  ylab(paste("PC2, ", pve_pc2, "% variance explained", sep =""))
pca_plot


head(test3)
head(test)
%>%
  mutate_at(vars(starts_with("lib")), ~(scale(.) %>% as.vector))

head(test2)


... or you could just do dat[columns] <- scale(dat[columns]), which has worked consistently for the past 20 years ;-) –

dat2 <- dat %>% mutate_at(c("y", "z"), ~(scale(.) %>% as.vector))
dat2
test2 = test[, -c(1,2,3)]

test2 = as.matrix(test2)

scale(test2)

%>% mutate_at(vars(starts_with("lib")), funs(c(scale(.))))

head(test2)
     mutate_at(c(3,6), funs(c(scale(.))))



frags %>% ggplot(., aes(x = start, y = ratio))

head(frag_count)

frags %>% pivot_wider(names_from = library_id)
test2

test2$current_dx
libraries_full$institution

  names(libraries_full)
ls()
head(test)
  group_by

  pivot_wider(names_from = station, values_from = seen)




head(frag_count)
