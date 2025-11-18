# Le génome est constitué de n paires de chromosomes sous forme d'une matrice
# n lignes par 2 colonnes
# les allèles de la composante paternelle se trouvent dans la 1ère colonne
# les allèles de la composante maternelle se trouvent dans la 2ème colonne

# generate a pool of chr_count chromatides from an initial pool that has to be defined
# method define if the alleles have to be choosen randomly or in a parental way among 
# the available pool of alleles
# - is.parental = FALSE : we only sample among the pool and sort in decreasing order
# - is.parental = TRUE  : for each chromosome we choose the paternal or maternal chromatide
# the pool in the parental method must correspond to maternal or paternal chromosomes 
generate_autosome_alleles <- function(pool_alleles, autosome_count) {

  autosome_alleles <- pool_alleles |> 
    sample(size = autosome_count, replace = FALSE) |> 
    sort(decreasing = TRUE)
  
  return(autosome_alleles)

}

# generate sex allele TRUE = X, FALSE = Y
generate_sex_allele <- function(is.female) {
   
  sex_allele <- sample(c(TRUE, FALSE), size = 1) | is.female

}

generate_male_sex_allele   <- function() {

  male_sex_allele <- generate_sex_allele(is.female = FALSE)

}

generate_female_sex_allele <- function() {

  female_sex_allele <- generate_sex_allele(is.female = TRUE)

}

generate_chromatides <- function(pool_alleles = 1:1e9, autosome_count = 22, is.female) {

  autosome_alleles <- generate_autosome_alleles(
    pool_alleles = pool_alleles,
    autosome_count = autosome_count
  )

  sex_allele <- generate_sex_allele(is.female)

  return(c(autosome_alleles, sex_allele))

}

generate_chromatides_from_male <- function() {

  chromatides_from_male <- generate_chromatides(is.female = FALSE)
  
}

generate_chromatides_from_female <- function() {

  chromatides_from_female <- generate_chromatides(is.female = TRUE)
  
}

generate_genome <- function(random = FALSE, is.daughter = TRUE) {

  chromatides <- c(generate_chromatides_from_male(), generate_chromatides_from_female()) 
  genome <- matrix(chromatides, ncol = 2, byrow = FALSE)
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
  
}

generate_genome_garcon <- function() {
  
  genome_garcon <- generate_genome(is.daughter = FALSE)

}

is.male <- function(genome) {

  n <- nrow(genome)
  return(xor(genome[n, 1], genome[n, 2]))

}


# la fonction new_alleles génère les allèles issus du père, ou de la mère.
# la fonction new_genome rassemble les 2 composantes issues de new_alleles

new_chromatides <- function(genome, set.sex = FALSE, male = TRUE) {

  n <- nrow(genome)
  chromatides_chosen <- sample(x = 1:2, size = n, replace = TRUE, prob = c(0.5, 0.5))
  chromatides_id <- sapply(1:n, \(i) { 2 * (i-1) + chromatides_chosen[i] })
  # I'm using transpose function to convert matrix into vector by rows
  new_chromatides <- matrix(t(genome)[chromatides_id], ncol = 1)
  if ( set.sex & is.male(genome)) {
    if (male) {
      new_chromatides[n] <- 0
    } else {
      new_chromatides[n] <- 1
    }
  }

  return(new_chromatides)

}

new_genome <- function(genome_father, genome_mother) {
    
  new_chromatides_from_father <- new_chromatides(genome_father)
  new_chromatides_from_mother <- new_chromatides(genome_mother)
  new_genome <- cbind( new_chromatides_from_father, new_chromatides_from_mother )

}

new_male_genome <- function(genome_father, genome_mother) {

  new_chromatides_from_father <- new_chromatides(genome_father, set.sex = TRUE, male = TRUE)
  new_chromatides_from_mother <- new_chromatides(genome_mother)
  new_genome <- cbind( new_chromatides_from_father, new_chromatides_from_mother )

}

new_female_genome <- function(genome_father, genome_mother) {

  new_chromatides_from_father <- new_chromatides(genome_father, set.sex = TRUE, male = FALSE)
  new_chromatides_from_mother <- new_chromatides(genome_mother)
  new_genome <- cbind( new_chromatides_from_father, new_chromatides_from_mother )

}

# proportions d'allèles du parent dans la composante associée de l'enfant
# male = TRUE : on compare les alleles issus du père
# male = TRUE : on compare les alleles issus de la mère
percentage_identical_alleles <- function(genome, other_genome, male = TRUE) {
 
  if (male) {
    col <- 1
  } else {
    col <- 2
  }
    
  alleles_to_compare <- as.vector(genome[, col])
  alleles_other_genome <- as.vector(other_genome)
  genome_size <- length(alleles_other_genome)
  percentage <- sum(alleles_to_compare %in% alleles_other_genome) / genome_size

}

simulation <- function(count, size) {

  list_genomes <- vector()
  for (i in 1:count) {
        
    prop.oncle <- numeric(size)
    for (j in 1:size) {

      # je vais simuler la naissance de 2 fils issus des deux grands-parents gp1 et gm1
      # donc possiblement 2 pères

      genome_p1_1 <- new_male_genome(genome_gp1, genome_gm1)
      genome_p1_2 <- new_male_genome(genome_gp1, genome_gm1)

      # je simule la naissance d'une fille issue des deux grands-parents gp2 et gm2
      # possiblement une mère

      genome_m2 <- new_female_genome(genome_gp2, genome_gm2)

      # je simule la naissance d'un individu issu de gp1_1 et gm2

      genome_p1_1_m2 <- new_male_genome(genome_p1_1, genome_m2)

      # proportions d'allèles de p1_1, p1_2 et m2 dans p1_1_m2
      # prop.pere[j] <- percentage_identical_alleles(genome_p1_1_m2, genome_p1_1)
      # prop.mere[j] <- percentage_identical_alleles(genome_p1_1_m2, genome_m2, male = FALSE)
      prop.oncle[j] <- percentage_identical_alleles(genome_p1_1_m2, genome_p1_2)
      if (prop.oncle[j] == 0.5) {
        list_genomes <- rbind(list_genomes, list(genome_p1_1, genome_p1_2, genome_m2, genome_p1_1_m2))
      }
    }
    saveRDS(prop.oncle, file = paste0("data/df", i, ".RDS"))
    }
  tib_genomes <- tibble::as.tibble(list_genomes)
  names(tib_genomes) <- c("genome.pere", "genome.oncle", "genome.mere", "genome.enfant")
  saveRDS(tib_genomes, file = paste0("tib_genomes.RDS"))

}

