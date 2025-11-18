source("R/functions.R")

# Génomes grands-parents paternels (suffixés 1) et maternels (suffixés 2)

genome_gp1 <- generate_genome_garcon()
genome_gm1 <- generate_genome_fille()
genome_gp2 <- generate_genome_garcon()
genome_gm2 <- generate_genome_fille()

simulation(count = 1000, size = 100000)

rds_files <- list.files(path = "data", pattern = "*.RDS" , full.names = TRUE)
prop.oncle = readRDS(rds_files[1])
df = data.frame(prop.oncle)
count <- length(rds_files)
for (i in 2:count) {

    prop.oncle <- readRDS(rds_files[i])
    df <- rbind(df, data.frame(prop.oncle = prop.oncle))

}

hist(df$prop.oncle, breaks = 12)
summary(df$prop.oncle)
#library(ggplot2)
#df |>
#  ggplot( aes(x=prop.oncle)) +
  #  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)

