# migration modified to match python

#dominance or overdominance
fitness_mode <- function(fitness, pop){
  if(fitness == 1){ #dominance
    locus <- matrix(1, nrow = 2, ncol = pop)
  }
  if(fitness == 2){ #overdominance
    locus <- matrix(1, nrow = 2, ncol = pop)
    locus[2,] <- rep(0, times = pop)
  }
  return(locus)
}

#only two allele can exist, either 1 or 0
mutation_bi <- function(locus, probability){
  x <- matrix(sample(1:1/probability, nrow(locus)*ncol(locus), replace = T), nrow = nrow(locus), ncol = ncol(locus))
  for(i in 1:nrow(locus)){
    for(j in 1:ncol(locus)){
      if(x[i, j] > 1){
        locus[i, j] <- locus[i, j]
      }else{
        locus[i, j] <- abs(locus[i, j] - 1)
      }
    }
  }
  return(locus)
}

#will random mutate
mutation_pro <- function(locus, probability){
  x <- matrix(sample(1:1/probability, nrow(locus)*ncol(locus), replace = T), nrow = nrow(locus), ncol = ncol(locus))
  for(i in 1:nrow(locus)){
    for(j in 1:ncol(locus)){
      if(x[i, j] > 1){
        locus[i, j] <- locus[i, j]
      }else{
        locus[i, j] <- sample(1:50, 1)/100
      }
    }
  }
  return(locus)
}

#can or can not encounter mates
cross <- function(encounter, sex){
  mate <- vector("numeric", length(sex))
  for(i in 1:length(sex)){
    if(sex[i] == 0){
      n <- ifelse(encounter < length(sex) - 1, encounter, length(sex) - 1)
      indiv <- sample(1:(length(sex) - 1), n, replace = T)
      indiv[indiv >= i] <- indiv[indiv >= i] + 1
      while(n > 0){
        if(sex[indiv[n]] == 1){
          mate[i] <- indiv[n]
          break
        }
        n <- n - 1
      }
    }
  }
  return(mate)
}

#which reproduction method female use
repro_method <- function(sex, locus2, mate){
  sexual <- vector("numeric", 0)
  par_sex <- vector("numeric", 0)
  par <- vector("numeric", 0)
  for(i in 1:length(sex)){
    if(sex[i] == 0 && mate[i] != 0){
      if(sample(1:100, 1)/100 > sum(locus2[,i])/2){
        sexual <- c(sexual, i)
      }else{
        par_sex <- c(par_sex, i)
      }
    }
    if((sex[i] == 0) && (mate[i] == 0) && (sample(1:100, 1)/100 <= sum(locus2[,i])/2)){
      par <- c(par, i)
    }
  }
  return(list(sexual = sexual, par_sex = par_sex, par = par))
}

#forming gamate
gamate <- function(locus1, locus2, locus3, recom12, recom23){
  x <- sample(0:1, 1)
  y <- sample(0:1/recom12, 1)
  z <- sample(0:1/recom23, 1)
  offspring_locus1 <- locus1[x + 1]
  if(y > 1){
    offspring_locus2 <- locus2[x + 1]
  }else{
    offspring_locus2 <- locus2[abs(x - 1) + 1]
  }
  if(z > 1){
    offspring_locus3 <- locus3[x + 1]
  }else{
    offspring_locus3 <- locus3[abs(x - 1) + 1]
  }
  return(c(offspring_locus1, offspring_locus2, offspring_locus3))
}

#selection based on fitness mode
selection <- function(fitness, locus_1, locus_2){
  live <- 0
  if(fitness == 1){
    if(locus_1 + locus_2 != 0){
      live <- 1
    }
  }else{
    if(locus_1 + locus_2 == 1){
      live <- 1
    }
  }
  return(live)
}

#reproduction
reproduction <- function(sexual, par_sex, par, Nsexual, Npar_sex, Npar, pop, locus1, locus2, locus3, recom12, recom23, mate, fitness){
  total <- length(sexual)*Nsexual + length(par_sex)*Npar_sex + length(par)*Npar
  size <- pop
  offlocus1 <- matrix(data = NA, nrow = 2, ncol = 0)
  offlocus2 <- matrix(data = NA, nrow = 2, ncol = 0)
  offlocus3 <- matrix(data = NA, nrow = 2, ncol = 0)
  offsex <- vector("numeric", 0)
  par_num <- 0
  if(total < pop){
    size <- total
  }
  i <- 1
  while(i <= size){
    type <- sample(1:total, 1)
    sex <- 100
    if(type <= length(sexual)*Nsexual){
      mother <- sexual[sample(1:length(sexual), 1)]
      mgamate <- gamate(locus1[,mother], locus2[,mother], locus3[,mother], recom12, recom23)
      fgamate <- gamate(locus1[,mate[mother]], locus2[,mate[mother]], locus3[,mate[mother]], recom12, recom23)
    }else if(type <= length(sexual)*Nsexual + length(par_sex)*Npar_sex){
      mother <- par_sex[sample(1:length(par_sex), 1)]
      mgamate <- gamate(locus1[,mother], locus2[,mother], locus3[,mother], recom12, recom23)
      fgamate <- gamate(locus1[,mate[mother]], locus2[,mate[mother]], locus3[,mate[mother]], recom12, recom23)
    }else{
      mother <- par[sample(1:length(par), 1)]
      mgamate <- gamate(locus1[,mother], locus2[,mother], locus3[,mother], recom12, recom23)
      fgamate <- mgamate
      sex <- 0
    }
    if(selection(fitness, mgamate[1], fgamate[1]) == 1){
      offlocus1 <- cbind(offlocus1, c(mgamate[1], fgamate[1]))
      offlocus2 <- cbind(offlocus2, c(mgamate[2], fgamate[2]))
      offlocus3 <- cbind(offlocus3, c(mgamate[3], fgamate[3]))
      if(sex == 0){
        offsex <- c(offsex, 0)
        par_num <- par_num + 1
      }else{
        offsex <- c(offsex, sample(0:1, 1))
      }
    }else{
      i <- i - 1
    }
    i <- i + 1
  }
  return(list(offlocus1, offlocus2, offlocus3, offsex, par_num))
}

args <- commandArgs(trailingOnly = TRUE)

# modified order to work with python parameters directly
mainpop <- as.numeric(args[1])
maxsubpop <- as.numeric(args[2])
migration_input <- as.numeric(args[3])
migration <- 1 / migration_input
Rmut_input <- as.numeric(args[4])
Rmut <- 1 / Rmut_input
recom12 <- as.numeric(args[5])
recom23 <- as.numeric(args[6])
main_encounter <- as.numeric(args[7])
sub_encounter <- max(1, min(floor(0.04 * maxsubpop), main_encounter %/% 2))
Nsexual <- as.numeric(args[8])
fitness <- (as.numeric(args[9]) + 1)
generation <- as.numeric(args[10])
Npar_input <- as.numeric(args[11])
Npar <- as.integer(Nsexual * Npar_input + 0.5)
Npar_sex_input <- as.numeric(args[12])
Npar_sex <- as.integer(Nsexual * (1 - Npar_sex_input) + 0.5)
simulation_id <- args[13]
run <- 96

run_simulation <- function(run_id) {
  set.seed(sample.int(1e6, 1) + run_id)

  main_sex <- sample(0:1, mainpop, replace = TRUE)
  sub_sex <- vector("numeric", 0)
  subpop <- 0
  main_locus1 <- fitness_mode(fitness, mainpop)
  main_locus2 <- matrix(data = 0, nrow = 2, ncol = mainpop)
  main_locus3 <- matrix(data = 0, nrow = 2, ncol = mainpop)
  sub_locus1 <- matrix(data = NA, nrow = 2, ncol = 0)
  sub_locus2 <- matrix(data = NA, nrow = 2, ncol = 0)
  sub_locus3 <- matrix(data = NA, nrow = 2, ncol = 0)
  f_locus2_main <- vector("numeric", generation)
  f_locus3_main <- vector("numeric", generation)
  f_locus2_sub <- vector("numeric", generation)
  f_locus3_sub <- vector("numeric", generation)
  par_num_main <- vector("numeric", generation)
  par_num_sub <- vector("numeric", generation)
  popsize_sub <- vector("numeric", generation)
  sexratio_main <- vector("numeric", generation)
  sexratio_sub <- vector("numeric", generation)

  for(g in 1:generation){
    par_sub <- 0
    main_locus1 <- mutation_bi(main_locus1, Rmut)
    main_locus2 <- mutation_bi(main_locus2, Rmut)
    main_locus3 <- mutation_bi(main_locus3, Rmut)
    main_mate <- cross(main_encounter, main_sex)
    main_repro <- repro_method(main_sex, main_locus2, main_mate)
    if(fitness == 2 && length(main_repro$sexual) + length(main_repro$par_sex) == 0){
      main_offspring <- list()
    }else{
      main_offspring <- reproduction(main_repro$sexual, main_repro$par_sex, main_repro$par, Nsexual, Npar_sex, Npar, mainpop, main_locus1, main_locus2, main_locus3, recom12, recom23, main_mate, fitness)
    }
    if(subpop != 0){
      sub_locus1 <- mutation_bi(sub_locus1, Rmut)
      sub_locus2 <- mutation_bi(sub_locus2, Rmut)
      sub_locus3 <- mutation_bi(sub_locus3, Rmut)
      sub_mate <- cross(sub_encounter, sub_sex)
      sub_repro <- repro_method(sub_sex, sub_locus2, sub_mate)
      if(fitness == 2 && (length(sub_repro$sexual) + length(sub_repro$par_sex) == 0)){
        sub_offspring <- list()
      }else{
        sub_offspring <- reproduction(sub_repro$sexual, sub_repro$par_sex, sub_repro$par, Nsexual, Npar_sex, Npar, maxsubpop, sub_locus1, sub_locus2, sub_locus3, recom12, recom23, sub_mate, fitness)
        par_sub <- sub_offspring[[5]]
      }
    }
    main_sex <- main_offspring[[4]]
    main_locus1 <- main_offspring[[1]]
    main_locus2 <- main_offspring[[2]]
    main_locus3 <- main_offspring[[3]]
    if(subpop != 0 && length(sub_offspring) != 0){
      sub_sex <- sub_offspring[[4]]
      sub_locus1 <- sub_offspring[[1]]
      sub_locus2 <- sub_offspring[[2]]
      sub_locus3 <- sub_offspring[[3]]
    }
    subpop <- length(sub_sex)
    mig_main <- vector("numeric", 0)
    for(i in 1:length(main_sex)){
      if(sample(1:1/migration, 1) <= 1){
        sub_sex <- c(sub_sex, main_sex[i])
        sub_locus1 <- cbind(sub_locus1, main_locus1[,i])
        sub_locus2 <- cbind(sub_locus2, main_locus2[,i])
        sub_locus3 <- cbind(sub_locus3, main_locus3[,i])
        mig_main <- cbind(mig_main, i)
      }
    }
    mig_sub <- vector("numeric", 0)
    if(subpop != 0){
      for(i in 1:subpop){
        if(sample(1:1/migration, 1) <= 1){
          main_sex <- c(main_sex, sub_sex[i])
          main_locus1 <- cbind(main_locus1, sub_locus1[,i])
          main_locus2 <- cbind(main_locus2, sub_locus2[,i])
          main_locus3 <- cbind(main_locus3, sub_locus3[,i])
          mig_sub <- cbind(mig_sub, i)
        }
      }
    }
    if(length(mig_main) != 0){
      main_sex <- main_sex[-mig_main]
      main_locus1 <- main_locus1[,-mig_main]
      main_locus2 <- main_locus2[,-mig_main]
      main_locus3 <- main_locus3[,-mig_main]
    }
    if(length(mig_sub) != 0){
      sub_sex <- sub_sex[-mig_sub]
      sub_locus1 <- sub_locus1[,-mig_sub]
      sub_locus2 <- sub_locus2[,-mig_sub]
      sub_locus3 <- sub_locus3[,-mig_sub]
    }
    if(length(main_sex) > mainpop){
      main_sex <- main_sex[-(mainpop + 1:length(main_sex))]
      main_locus1 <- main_locus1[,-(mainpop + 1:length(main_locus1))]
      main_locus2 <- main_locus2[,-(mainpop + 1:length(main_locus2))]
      main_locus3 <- main_locus3[,-(mainpop + 1:length(main_locus3))]
    }
    if(length(sub_sex) > maxsubpop){
      sub_sex <- sub_sex[-(maxsubpop + 1:length(sub_sex))]
      sub_locus1 <- sub_locus1[,-(maxsubpop + 1:length(sub_locus1))]
      sub_locus2 <- sub_locus2[,-(maxsubpop + 1:length(sub_locus2))]
      sub_locus3 <- sub_locus3[,-(maxsubpop + 1:length(sub_locus3))]
    }
    subpop <- length(sub_sex)
    f_locus2_main[g] <- sum(main_locus2)/ncol(main_locus2)/2
    f_locus3_main[g] <- sum(main_locus3)/ncol(main_locus3)/2
    par_num_main[g] <- main_offspring[[5]]
    popsize_sub[g] <- subpop
    sexratio_main[g] <- 1 - (sum(main_sex)/length(main_sex))
    if(subpop != 0){
      f_locus2_sub[g] <- sum(sub_locus2)/ncol(sub_locus2)/2
      f_locus3_sub[g] <- sum(sub_locus3)/ncol(sub_locus3)/2
      par_num_sub[g] <- par_sub
      sexratio_sub[g] <- 1 - (sum(sub_sex)/length(sub_sex))
    }else{
      f_locus2_sub[g] <- NA
      f_locus3_sub[g] <- NA
      par_num_sub[g] <- NA
      sexratio_sub[g] <- NA
    }
  }

  return(list(f_locus2_main=f_locus2_main, f_locus2_sub=f_locus2_sub, f_locus3_main=f_locus3_main, f_locus3_sub=f_locus3_sub, par_num_main=par_num_main, par_num_sub=par_num_sub, popsize_sub=popsize_sub, sexratio_main=sexratio_main, sexratio_sub=sexratio_sub))

}

library(parallel)
library(parallelly)
results <- mclapply(1:run, run_simulation, mc.cores = availableCores())

all_locus2_main <- do.call(rbind, lapply(results, function(x) x$f_locus2_main))
all_locus2_sub <- do.call(rbind, lapply(results, function(x) x$f_locus2_sub))
all_locus3_main <- do.call(rbind, lapply(results, function(x) x$f_locus3_main))
all_locus3_sub <- do.call(rbind, lapply(results, function(x) x$f_locus3_sub))
all_par_num_main <- do.call(rbind, lapply(results, function(x) x$par_num_main))
all_par_num_sub <- do.call(rbind, lapply(results, function(x) x$par_num_sub))
all_popsize_sub <- do.call(rbind, lapply(results, function(x) x$popsize_sub))
all_sexratio_main <- do.call(rbind, lapply(results, function(x) x$sexratio_main))
all_sexratio_sub <- do.call(rbind, lapply(results, function(x) x$sexratio_sub))

# results
output_dir <- file.path("/gpfs/home/blukacsy/true/graphs_R", paste0("sim", simulation_id))
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# plotting

# locus 2
png(filename = file.path(output_dir, "locus2.png"), width = 1200, height = 600, type = "cairo")

par(mfrow = c(1, 2))
plot(x = c(0,generation), y = c(0,1), type = "n", xlab = "generation", ylab = "frequency", main = "mainpopulation parthenogenesis locus")
for(i in 1:run){
  lines(all_locus2_main[i,], col = i)
}

plot(x = c(0,generation), y = c(0,1), type = "n", xlab = "generation", ylab = "frequency", main = "subpopulation parthenogenesis locus")
for(i in 1:run){
  lines(all_locus2_sub[i,], col = i)
}

dev.off()

# locus 3
png(filename = file.path(output_dir, "locus3.png"), width = 1200, height = 600, type = "cairo")

par(mfrow = c(1, 2))
plot(x = c(0,generation), y = c(0,1), type = "n", xlab = "generation", ylab = "frequency", main = "mainpopulation neutral locus")
for(i in 1:run){
  lines(all_locus3_main[i,], col = i)
}

plot(x = c(0,generation), y = c(0,1), type = "n", xlab = "generation", ylab = "frequency", main = "subpopulation neutral locus")
for(i in 1:run){
  lines(all_locus3_sub[i,], col = i)
}

dev.off()

# subpop count
png(filename = file.path(output_dir, "subpop_size.png"), width = 800, height = 600, type = "cairo")

plot(x = c(0,generation), y = c(0,maxsubpop), type = "n", xlab = "generation", ylab = "Number of flies", main = "Subpopulation Size")
for(i in 1:run){
  lines(all_popsize_sub[i,], col = i)
}

dev.off()

# tables
write.table(all_locus2_main, file = file.path(output_dir, "all_locus2_main.txt"), sep = ",", row.names = FALSE)
write.table(all_locus2_sub,  file = file.path(output_dir, "all_locus2_sub.txt"),  sep = ",", row.names = FALSE)
write.table(all_locus3_main, file = file.path(output_dir, "all_locus3_main.txt"), sep = ",", row.names = FALSE)
write.table(all_locus3_sub,  file = file.path(output_dir, "all_locus3_sub.txt"),  sep = ",", row.names = FALSE)
write.table(all_par_num_main, file = file.path(output_dir, "all_par_num_main.txt"), sep = ",", row.names = FALSE)
write.table(all_par_num_sub,  file = file.path(output_dir, "all_par_num_sub.txt"),  sep = ",", row.names = FALSE)
write.table(all_popsize_sub,  file = file.path(output_dir, "all_popsize_sub.txt"),  sep = ",", row.names = FALSE)
write.table(all_sexratio_main, file = file.path(output_dir, "all_sexratio_main.txt"), sep = ",", row.names = FALSE)
write.table(all_sexratio_sub,  file = file.path(output_dir, "all_sexratio_sub.txt"),  sep = ",", row.names = FALSE)

# summary
parthe_main <- vector("numeric", 0)
for(i in 1:nrow(all_locus2_main)){
  if(1 %in% all_locus2_main[i,]){
    parthe_main <- c(parthe_main, i)
  }
}

parthe_sub <- vector("numeric", 0)
for(i in 1:nrow(all_locus2_sub)){
  for(j in 1:ncol(all_locus2_sub)){
    if(is.na(all_locus2_sub[i,j])){
      next
    }
    if(all_locus2_sub[i,j] == 1 && all_popsize_sub[i, j] > 10){
      parthe_sub <- c(parthe_sub, i)
      break
    }
  }
}

neutral_main <- vector("numeric", 0)
for(i in 1:nrow(all_locus3_main)){
  if(1 %in% all_locus3_main[i,]){
    neutral_main <- c(neutral_main, i)
  }
}

neutral_sub <- vector("numeric", 0)
for(i in 1:nrow(all_locus3_sub)){
  for(j in 1:ncol(all_locus3_sub)){
    if(is.na(all_locus3_sub[i,j])){
      next
    }
    if(all_locus3_sub[i,j] == 1 && all_popsize_sub[i, j] > 10){
      neutral_sub <- c(neutral_sub, i)
      break
    }
  }
}

nosub <- vector("numeric", 0)
for(i in 1:nrow(all_popsize_sub)){
  if(!(maxsubpop %in% all_popsize_sub[i,])){
    nosub <- c(nosub, i)
  }
}

interact <- vector("numeric", 0)
for(i in 1:nrow(all_popsize_sub)){
  if(maxsubpop %in% all_par_num_sub[i,]){
    interact <- c(interact, i)
    next
  }
  for(j in 1:ncol(all_popsize_sub)){
    if(!is.na(all_par_num_sub[i,j]) && all_par_num_sub[i,j] != 0){
      expmale <- all_popsize_sub[i,j] - all_par_num_sub[i,j]
      if((round((1 - all_sexratio_sub[i,j]) * all_popsize_sub[i,j]) == expmale) && expmale != 0){
        interact <- c(interact, i)
        break
      }
    }
  }
}

max_len <- max(length(parthe_main), length(parthe_sub), length(neutral_main), length(neutral_sub), length(nosub), length(interact))

pad_vector <- function(vec, max_len) {
  c(vec, rep(NA, max_len - length(vec)))
}

parthe_main_pad <- pad_vector(parthe_main, max_len)
parthe_sub_pad <- pad_vector(parthe_sub, max_len)
neutral_main_pad <- pad_vector(neutral_main, max_len)
neutral_sub_pad <- pad_vector(neutral_sub, max_len)
nosub_pad <- pad_vector(nosub, max_len)
interact_pad <- pad_vector(interact, max_len)

summary_data <- data.frame(
  "parth_main_fixed_run" = parthe_main_pad,
  "parth_sub_fixed_run" = parthe_sub_pad,
  "neutral_main_fixed_run" = neutral_main_pad,
  "neutral_sub_fixed_run" = neutral_sub_pad,
  "no_max_sub_run" = nosub_pad,
  "interact_run" = interact_pad
)

write.table(summary_data, file = file.path(output_dir, "summary.txt"), sep = ",", row.names = FALSE, col.names = TRUE)
