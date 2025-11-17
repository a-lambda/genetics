# Le génome est constitué de n paires de chromosomes sous forme d'une matrice
# n lignes par 2 colonnes
# les allèles de la composante paternelle se trouvent dans la 1ère colonne
# les allèles de la composante maternelle se trouvent dans la 2ème colonne

# Génération d'un génome


CHR_COUNT <- 22
POOL_ALLELES <- 1:1e9

choose_autosomes_alleles <- function(pool_alleles = POOL_ALLELES, chr_count = CHR_COUNT) {
    
    autosomes_alleles <- pool_alleles |> 
        sample(size = chr_count, replace = FALSE) |> 
        sort(decreasing = TRUE)
    
    return(autosomes_alleles)

}

choose_sex_allele <- function(is.female) {
   
    sex_allele <- sample(c(TRUE, FALSE), size = 1) | is.female
    return(sex_allele)

}

choose_male_sex_allele   <- function() { choose_sex_allele(is.female = FALSE) }
choose_female_sex_allele <- function() { choose_sex_allele(is.female = TRUE) }

generate_chromatides <- function(pool_alleles = POOL_ALLELES, chr_count = CHR_COUNT, is.female) {

    autosomes_alleles <- choose_autosomes_alleles(
        pool_alleles = pool_alleles,
        chr_count = chr_count
    )
    sex_allele <- choose_sex_allele(is.female)

    return(c(autosomes_alleles, sex_allele))

}

generate_chromatides_from_male <- function() {

    chromatides_from_male <- generate_chromatides(is.female = FALSE)
    
    return(chromatides_from_male)
  
}

generate_chromatides_from_female <- function() {

    chromatides_from_female <- generate_chromatides(is.female = TRUE)
  
    return(chromatides_from_female)
   
}

generate_genome <- function(random = FALSE, is.daughter = TRUE) {

    genome <- matrix(
        c(generate_chromatides_from_male(), generate_chromatides_from_female()),
        ncol = 2,
        byrow = FALSE
    )
    
    if (random) { 
      return(genome) 
    } else {
      chr_count <- nrow(genome)
      genome[chr_count, 1] <- as.numeric(is.daughter)
      return(genome)
    }
}

generate_genome_fille <- function() {
  
    genome_fille <- generate_genome()
    return(genome_fille)
  
}

generate_genome_garcon <- function() {
  
    genome_garcon <- generate_genome(is.daughter = FALSE)
    return(genome_garcon)

}

is.male <- function(genome, chr_count = CHR_COUNT) {

    return(xor(genome[chr_count + 1, 1], genome[chr_count + 1, 2]))

}

# la fonction new_alleles génère les allèles issus du père, ou de la mère.
# la fonction new_genome rassemble les 2 composantes issues de new_alleles

new_alleles <- function(genome) {

    nr <- nrow(genome)
    alleles_chosen <- sample(x = 1:2, size = nr, replace = T, prob = c(0.5, 0.5))
    alleles_id <- sapply(1:nr, \(i) {2 * (i-1) + alleles_chosen[i]})
    # I'm using transpose function to convert matrix into vector by rows
    new_alleles <- matrix(t(genome)[alleles_id], ncol = 1)
    return (new_alleles)

}

new_genome <- function(genome_father, genome_mother) {
    
    new_alleles_from_father <- new_alleles(genome_father)
    new_alleles_from_mother <- new_alleles(genome_mother)
    new_genome <- cbind(
        new_alleles_from_father,
        new_alleles_from_mother
        )

    return(new_genome)

}

# proportions d'allèles du parent dans la composante associée de l'enfant

percentage_of_alleles_parent <- function(genome_child, genome_parent, male = TRUE) {
 
    if (male) {
        col <- 1
    } else {
        col <- 2
    }
    
    alleles_child_to_compare <- as.vector(genome_child[, col])
    alleles_parent <- as.vector(genome_parent)
    genome_size <- length(alleles_parent)
    percentage <- sum(alleles_child_to_compare %in% alleles_parent) / genome_size

    return(percentage)

}

simulation <- function(count, size) {

    list_genomes <- vector()
    for (i in 1:count) {
        
        prop.oncle <- numeric(size)
        for (j in 1:size) {

# je vais simuler la naissance de 2 fils issus des deux grands-parents gp1 et gm1
# pas de crossover
# donc possiblement 2 pères

            genome_p1_1 <- new_genome(genome_gp1, genome_gm1)
            genome_p1_2 <- new_genome(genome_gp1, genome_gm1)

# je simule la naissance d'une fille issue des deux grands-parents gp2 et gm2
# possiblement une mère

            genome_m2 <- new_genome(genome_gp2, genome_gm2)

# je simule la naissance d'un individu issu de p1_1 et m2

            genome_p1_1_m2 <- new_genome(genome_p1_1, genome_m2)

# genomes obtenus

            # genome_p1_1
            # genome_p1_2
            # genome_m2
            # genome_p1_1_m2


# proportions d'allèles de p1_1, p1_2 et m2 dans p1_1_m2

            # prop.pere[j] <- percentage_of_alleles_parent(genome_p1_1_m2, genome_p1_1)
            prop.oncle[j] <- percentage_of_alleles_parent(genome_p1_1_m2, genome_p1_2)
            if (prop.oncle[j] == 0.5) {
                list_genomes <- rbind(list_genomes, list(genome_p1_1, genome_p1_2, genome_m2, genome_p1_1_m2))
            }
            # prop.mere[j] <- percentage_of_alleles_parent(genome_p1_1_m2, genome_m2, male = FALSE)

        }
        saveRDS(prop.oncle, file = paste0("data/df", i, ".RDS"))
    }
    tib_genomes <- tibble::as.tibble(list_genomes)
    names(tib_genomes) <- c("genome.pere", "genome.oncle", "genome.mere", "genome.enfant")
    saveRDS(tib_genomes, file = paste0("tib_genomes.RDS"))

}

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
