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

#cin
mainpop <- 500
maxsubpop <- 50
fitness <- 1
Rmut <- 1/10000
main_encounter <- 4
sub_encounter <- min(0.04*maxsubpop, main_encounter/2)
Nsexual <- 10
Npar_sex <- round(Nsexual * 0.75)
Npar <- 2
recom12 <- 0.5
recom23 <- 0.5
migration <- 1/10/mainpop
generation <- 1000
run <- 50
all_locus2_main <- matrix(0, nrow = run, ncol = generation)
all_locus2_sub <- matrix(0, nrow = run, ncol = generation)
all_locus3_main <- matrix(0, nrow = run, ncol = generation)
all_locus3_sub <- matrix(0, nrow = run, ncol = generation)
all_par_num_main <- matrix(0, nrow = run, ncol = generation)
all_par_num_sub <- matrix(0, nrow = run, ncol = generation)
all_popsize_sub <- matrix(0, nrow = run, ncol = generation)
all_sexratio_main <- matrix(0, nrow = run, ncol = generation)
all_sexratio_sub <- matrix(0, nrow = run, ncol = generation)

t <- Sys.time()

for(n in 1:run){
  #initial state
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
    main_sex <- vector("numeric", 0)
    sub_sex <- vector("numeric", 0)
    main_locus1 <- matrix(data = NA, nrow = 2, ncol = 0)
    main_locus2 <- matrix(data = NA, nrow = 2, ncol = 0)
    main_locus3 <- matrix(data = NA, nrow = 2, ncol = 0)
    sub_locus1 <- matrix(data = NA, nrow = 2, ncol = 0)
    sub_locus2 <- matrix(data = NA, nrow = 2, ncol = 0)
    sub_locus3 <- matrix(data = NA, nrow = 2, ncol = 0)
    for(i in 1:length(main_offspring[[4]])){
      if(sample(1:1/migration, 1) > 1){
        main_sex <- c(main_sex, main_offspring[[4]][i])
        main_locus1 <- cbind(main_locus1, main_offspring[[1]][,i])
        main_locus2 <- cbind(main_locus2, main_offspring[[2]][,i])
        main_locus3 <- cbind(main_locus3, main_offspring[[3]][,i])
      }else{
        sub_sex <- c(sub_sex, main_offspring[[4]][i])
        sub_locus1 <- cbind(sub_locus1, main_offspring[[1]][,i])
        sub_locus2 <- cbind(sub_locus2, main_offspring[[2]][,i])
        sub_locus3 <- cbind(sub_locus3, main_offspring[[3]][,i])
      }
    }
    if(subpop != 0 && length(sub_offspring) != 0){
      i <- length(sub_offspring[[4]])
      while(i > 0){
        if(sample(1:1/migration, 1) > 1){
          sub_sex <- c(sub_sex, sub_offspring[[4]][i])
          sub_locus1 <- cbind(sub_locus1, sub_offspring[[1]][,i])
          sub_locus2 <- cbind(sub_locus2, sub_offspring[[2]][,i])
          sub_locus3 <- cbind(sub_locus3, sub_offspring[[3]][,i])
        }else{
          main_sex <- c(main_sex, sub_offspring[[4]][i])
          main_locus1 <- cbind(main_locus1, sub_offspring[[1]][,i])
          main_locus2 <- cbind(main_locus2, sub_offspring[[2]][,i])
          main_locus3 <- cbind(main_locus3, sub_offspring[[3]][,i])
        }
        i <- i - 1
      }
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
  all_locus2_main[n,] <- f_locus2_main
  all_locus2_sub[n,] <- f_locus2_sub
  all_locus3_main[n,] <- f_locus3_main
  all_locus3_sub[n,] <- f_locus3_sub
  all_par_num_main[n,] <- par_num_main
  all_par_num_sub[n,] <- par_num_sub
  all_popsize_sub[n,] <- popsize_sub
  all_sexratio_main[n,] <- sexratio_main
  all_sexratio_sub[n,] <- sexratio_sub
  cat("Run", n, "\n")
}

t <- Sys.time() - t

par(mfrow = c(1, 2))
plot(x = c(0,generation), y = c(0,1), type = "n", xlab = "generation", ylab = "frequency", main = "mainpopulation parthenogenesis locus")
for(i in 1:run){
  lines(all_locus2_main[i,], col = i)
}

plot(x = c(0,generation), y = c(0,1), type = "n", xlab = "generation", ylab = "frequency", main = "subpopulation parthenogenesis locus")
for(i in 1:run){
  lines(all_locus2_sub[i,], col = i)
}


par(mfrow = c(1, 2))
plot(x = c(0,generation), y = c(0,1), type = "n", xlab = "generation", ylab = "frequency", main = "mainpopulation neutral locus")
for(i in 1:run){
  lines(all_locus3_main[i,], col = i)
}

plot(x = c(0,generation), y = c(0,1), type = "n", xlab = "generation", ylab = "frequency", main = "subpopulation neutral locus")
for(i in 1:run){
  lines(all_locus3_sub[i,], col = i)
}


write.table(all_locus2_main, file = "C:/Users/user/Desktop/fitness_1/data_500_mig_1/all_locus2_main.txt", sep = ",", row.names = F)
write.table(all_locus2_sub, file = "C:/Users/user/Desktop/fitness_1/data_500_mig_1/all_locus2_sub.txt", sep = ",", row.names = F)
write.table(all_locus3_main, file = "C:/Users/user/Desktop/fitness_1/data_500_mig_1/all_locus3_main.txt", sep = ",", row.names = F)
write.table(all_locus3_sub, file = "C:/Users/user/Desktop/fitness_1/data_500_mig_1/all_locus3_sub.txt", sep = ",", row.names = F)
write.table(all_par_num_main, file = "C:/Users/user/Desktop/fitness_1/data_500_mig_1/all_par_num_main.txt", sep = ",", row.names = F)
write.table(all_par_num_sub, file = "C:/Users/user/Desktop/fitness_1/data_500_mig_1/all_par_num_sub.txt", sep = ",", row.names = F)
write.table(all_popsize_sub, file = "C:/Users/user/Desktop/fitness_1/data_500_mig_1/all_popsize_sub.txt", sep = ",", row.names = F)
write.table(all_sexratio_main, file = "C:/Users/user/Desktop/fitness_1/data_500_mig_1/all_sexratio_main.txt", sep = ",", row.names = F)
write.table(all_sexratio_sub, file = "C:/Users/user/Desktop/fitness_1/data_500_mig_1/all_sexratio_sub.txt", sep = ",", row.names = F)

