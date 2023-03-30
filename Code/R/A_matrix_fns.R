#load("Data/networks/Yang_a_matrices/adjacent_matrix_list.Rdata")
#adjacent_matrix_key <- read.csv("Data/networks/Yang_a_matrices/adjacent_matrix_key.csv")


assign_A_matrix<-function(motif.id,adjacent_matrix,new_spp_ls,minimal = F){
  #A_matrix <- matrix(NA, 4, 4)
  A_matrix <- matrix(NA,nrow = nrow(adjacent_matrix),ncol = ncol(adjacent_matrix))
  
  # assign the diagonal entries - intraspecific competiton
  diagonal.elements <- diag(adjacent_matrix)
  diagonal.elements[diagonal.elements==0] <- -0.1
  diag(A_matrix) <- diagonal.elements
  
  # assign the interspecific competition between basal species
  # how many basal species
  basal.species.n <- length(which(diagonal.elements == -1)) 
  # all entries are sampled from the uniform distribution [-0.5,0]
  A_matrix_sub <- matrix(runif(basal.species.n^2,min=-0.5,max=0),basal.species.n,basal.species.n) 
  # change the diagonal entries of basal species to -1
  diag(A_matrix_sub) <--1  
  A_matrix[1:basal.species.n,1:basal.species.n] <- A_matrix_sub
  
  # assign the consumer-resource interaction coefficient
  for(i in (basal.species.n+1):nrow(adjacent_matrix)){ # basal species should not be a consumer
    # how many consumer links
    possible_links_from_species_i_to_j <- adjacent_matrix[1:(i-1), i]
    N_possible_links_from_species_i_to_j <- length(which(possible_links_from_species_i_to_j==-1))
    if(N_possible_links_from_species_i_to_j == 1){
      possible_links_from_species_i_to_j[possible_links_from_species_i_to_j==-1] <- -0.5
    }else{
      consumer_interaction_value <- c(-0.4, rep(-0.1/(N_possible_links_from_species_i_to_j-1), N_possible_links_from_species_i_to_j-1))
      #consumer_interaction_value <- c(-0.4, rep(-0.2/(N_possible_links_from_species_i_to_j-1), N_possible_links_from_species_i_to_j-1))
      #consumer_interaction_value <- c(-0.4, runif(N_possible_links_from_species_i_to_j-1,-0.4,0))
      permute_value <- sample(consumer_interaction_value, length(consumer_interaction_value))
      possible_links_from_species_i_to_j[possible_links_from_species_i_to_j==-1] <- permute_value
    }
    A_matrix[1:(i-1), i] <- possible_links_from_species_i_to_j
  }
  
  # assign the resource-consumer interaction coefficiens
  # generate the energy conversion rate matrix first
  conversion_matrix <- A_matrix
  conversion_matrix[!is.na(conversion_matrix)] <- 1
  conversion_matrix[is.na(conversion_matrix)] <- -0.2
  
  if(isFALSE(minimal)){
  if(motif.id %in% c(7,10)){ # add lower influence in circumstances where omnivory
    #conversion_matrix[4,1] <- -0.02
    conversion_matrix[new_spp_ls[new_spp_ls$tlvl == 4,2],new_spp_ls[new_spp_ls$tlvl == 1,2]] <- -0.02
    
  }else if(motif.id %in% c(8,9)){
      #conversion_matrix[3,1] <- -0.02
      conversion_matrix[new_spp_ls[new_spp_ls$tlvl == 3,2],new_spp_ls[new_spp_ls$tlvl == 1,2]] <- -0.02
      
    }else if(motif.id %in% c(11,13)){
        #conversion_matrix[3,1] <- -0.02
        conversion_matrix[new_spp_ls[new_spp_ls$tlvl == 3,2],new_spp_ls[new_spp_ls$tlvl == 1,2]] <- -0.02
        #conversion_matrix[4,1] <- -0.02
        conversion_matrix[new_spp_ls[new_spp_ls$tlvl == 4,2],new_spp_ls[new_spp_ls$tlvl == 1,2]] <- -0.02
        
      }else if(motif.id == 12){
          #conversion_matrix[4,1] <- -0.02
          conversion_matrix[new_spp_ls[new_spp_ls$tlvl == 4,2],new_spp_ls[new_spp_ls$tlvl == 1,2]] <- -0.02
          #conversion_matrix[4,2] <- -0.02
          conversion_matrix[new_spp_ls[new_spp_ls$tlvl == 4,2],new_spp_ls[new_spp_ls$tlvl == 2,2]] <- -0.02
        
          }else if(motif.id == 6){
            #conversion_matrix[4,2] <- -0.02
            conversion_matrix[new_spp_ls[new_spp_ls$tlvl == 4,2],new_spp_ls[new_spp_ls$tlvl == 2,2]] <- -0.02
          
            }else{
            conversion_matrix <- conversion_matrix
          }
  }
  
  if(isTRUE(minimal)){
    if(motif.id %in% c(3,6)){ # add lower influence in circumstances where omnivory
      #conversion_matrix[4,1] <- -0.02
      conversion_matrix[new_spp_ls[new_spp_ls$tlvl == 4,2],new_spp_ls[new_spp_ls$tlvl == 1,2]] <- -0.02
      
    }else if(motif.id %in% c(4,5)){
      #conversion_matrix[3,1] <- -0.02
      conversion_matrix[new_spp_ls[new_spp_ls$tlvl == 3,2],new_spp_ls[new_spp_ls$tlvl == 1,2]] <- -0.02
      
    }else if(motif.id %in% c(7,8)){
      #conversion_matrix[3,1] <- -0.02
      conversion_matrix[new_spp_ls[new_spp_ls$tlvl == 3,2],new_spp_ls[new_spp_ls$tlvl == 1,2]] <- -0.02
      #conversion_matrix[4,1] <- -0.02
      conversion_matrix[new_spp_ls[new_spp_ls$tlvl == 4,2],new_spp_ls[new_spp_ls$tlvl == 1,2]] <- -0.02
      
    }else if(motif.id == 10){
      #conversion_matrix[4,1] <- -0.02
      conversion_matrix[new_spp_ls[new_spp_ls$tlvl == 3,2],new_spp_ls[new_spp_ls$tlvl == 1,2]] <- -0.02
      
    }else{
      conversion_matrix <- conversion_matrix
    }
  }
  
  # to obtain resource-consumer interaction
  for(m in 1:ncol(adjacent_matrix)){
    for(n in 1:ncol(adjacent_matrix)){
      A_matrix[m,n] <- ifelse(is.na(A_matrix[m,n]),
                              conversion_matrix[m,n]*A_matrix[n,m],A_matrix[m,n])
    }
  }
  return(A_matrix)
}

assign_R_small<-function(adjacent_matrix,new_spp_ls){
  R <- rep(0,dim(new_spp_ls)[1]) #redefine R vector
  #adjacent_matrix<-adjacent_matrix_list[[motif.id]]
  R<-diag(adjacent_matrix)*(-1) #ensure diagonal positive
  basalN<-length(which(R==1))
  #R[1:(basalN)] <- sort(runif(1:basalN,0.1,1.0),decreasing = TRUE)
  
  R[(basalN+1):length(R)] <- -sort(runif((length(R)-basalN),0.001,0.1),decreasing = TRUE)
  #R[(basalN+1):length(R)] <- -sort(runif((length(R)-basalN),0,0.001),decreasing = TRUE) #sample growth rates between 0 and 0.001
  
   return(R)
}

assign_R_large<-function(adjacent_matrix,new_spp_ls){
  R <- rep(0,dim(new_spp_ls)[1]) #redefine R vector
  #adjacent_matrix<-adjacent_matrix_list[[motif.id]]
  R<-diag(adjacent_matrix)*(-1) #ensure diagonal positive
  basalN<-length(which(R==1))
  #R[1:(basalN)] <- sort(runif(1:basalN,0.1,1.0),decreasing = TRUE)
  
  #R[(basalN+1):length(R)] <- -sort(runif((length(R)-basalN),0.001,0.1),decreasing = TRUE)
  #R[(basalN+1):length(R)] <- -sort(runif((length(R)-basalN),0,0.001),decreasing = TRUE) #sample growth rates between 0 and 0.001
  R[(basalN+1):length(R)] <- -sort(runif((length(R)-basalN),0,0.01),decreasing = TRUE) #sample growth rates between 0 and 0.001
  
  return(R)
}

check_local_stability<-function(jacobian_matrix){
  eigen.values.real.part<-Re(eigen(jacobian_matrix)$values)
  local_stability<-ifelse(all(eigen.values.real.part<0),1,0)
  return(local_stability)
}

proportion_trophic <- function(n,first,second,third,fourth){
  
  if(sum(first,second,third,fourth) != 1.0){
    stop("Trophic level proportions do not sum to 1.0")
  }
  x <- seq_along(1:n)
  
  prop_first <- round(length(x)*first) #convert proportions to absolute count
  prop_second <- round(length(x)*second)
  prop_third <- round(length(x)*third)
  prop_fourth <- n - prop_first - prop_second - prop_third
  
  if(third == 0){ #for cases when no high trophic levels, ensure set to 0 and total species sum to n
    prop_second <- n - prop_first
    prop_third <- 0
  }
  
  if(fourth == 0){
    prop_third <- n - prop_first - prop_second
    prop_fourth <- 0
  }
  
  splt_ls <- split(x, rep(1:4, c(prop_first, 
                                 prop_second, 
                                 prop_third,
                                 prop_fourth))) #split number of species in to trophic levels
  out <- matrix(NA,nrow = n,ncol=2)
  for(i in 1:length(splt_ls)){
    out[min(splt_ls[[i]]):max(splt_ls[[i]]),] <- c(as.numeric(rep(names(splt_ls)[i],length(splt_ls[[i]]))),
                                                   splt_ls[[i]]) #recombine species in to matrix with their trophic levels labelled
  }
  return(as.data.frame(out) |> `colnames<-`(c("tlvl","species")))
}

generalise_A <- function(adjacent_matrix_ls,motif.id,key,new_spp_ls){
  
  ref_A <- adjacent_matrix_ls[[motif.id]] #extract reference interaction matrix
  ref_A[row(ref_A)==col(ref_A)] <- 0 #remove diagonal
  
  tmp <- data.frame(merge(key,expand.grid(key),by.x="species", by.y="species"),
                    "value"= c(ref_A)) |>
    `colnames<-`(c("species","target_tlvl","spp_tlvl","value")) #match trophic levels with their respective values in the reference key and reference matrix
  
  out_A <- matrix(NA,nrow=dim(new_spp_ls)[1],ncol =dim(new_spp_ls)[1]) #now compare the desired foodweb to the above reference
  for(i in 1:dim(new_spp_ls)[1]){
    for(j in 1:dim(new_spp_ls)[1]){
      if(i == j & subset(new_spp_ls,species == i)$tlvl == 1){
        out_A[i,j] = -1 #if a basal species, set high intraspecific competition
      }else{ #populate the matrix with the identified values once filtered to the correct interaction
        out_A[i,j] = subset(tmp,
                            spp_tlvl %in% subset(new_spp_ls,species == i)$tlvl
                            & target_tlvl %in% subset(new_spp_ls,species == j)$tlvl)$value[1]
      }
    }
  }
  return(as.matrix(out_A))
}

generate_foodwebs_small <- function(motif.id,adjacent_matrix,new_spp_ls,N,minimal = F){
  k<-0
  foodweb_out<-list()
  pb <-txtProgressBar(min = 0, max = N, style = 3)
  
  while(k<N){
    A_mat <- assign_A_matrix(motif.id,adjacent_matrix,new_spp_ls,minimal = minimal) # interaction coefficient matrix
    R_vec <- assign_R_small(adjacent_matrix,new_spp_ls)
    
    Neq <- tryCatch(base::solve(A_mat,-R_vec), 
                    error=function(err) rep(-1,dim(adjacent_matrix)[1])) #solve to equilbrium
    #Neq[Neq<0] <- Neq[Neq<0]*-0.01
    #Neq[Neq<0] <- Neq[Neq<0]*0+0.01
    J <- A_mat*Neq # Jacobian matrix
    L <- unlist(eigen(J)$values) # eigenvalues of the Jacobian matrix
    maxreL <- max(Re(L)) # max real part of Jacobian matrix
    ct <- all(Neq>=0) & maxreL < (-0.005) & maxreL > - 0.1 # check for local stability and feasibility
    #ct <- maxreL < (-0.005) & maxreL > - 0.1

    if(isTRUE(ct)){
      k<-k+1
      community<-list(J,R_vec,A_mat,Neq,maxreL,as.matrix(new_spp_ls))
      names(community)<-c("jacobian_matrix","R","A_matrix","Neq","maxreL","tlvl")
      foodweb_out[[k]]<-community   
      setTxtProgressBar(pb, k)
      
    }else{
      k<-k
    } 
  }    
  close(pb)
  
  return(foodweb_out)
}

generate_foodwebs_large <- function(motif.id,adjacent_matrix,new_spp_ls,N,minimal = F){
  k<-0
  foodweb_out<-list()
  pb <-txtProgressBar(min = 0, max = N, style = 3)
  while(k<N){
    A_mat <- assign_A_matrix(motif.id,adjacent_matrix,new_spp_ls,minimal = minimal) # interaction coefficient matrix
    R_vec <- assign_R_large(adjacent_matrix,new_spp_ls)
    
    Neq <- tryCatch(base::solve(A_mat,-R_vec), 
                    error=function(err) rep(-1,dim(adjacent_matrix)[1])) #solve to equilbrium
    #Neq[Neq<0] <- Neq[Neq<0]*-0.01
    #Neq[Neq<0] <- Neq[Neq<0]*0+0.01
    J <- A_mat*Neq # Jacobian matrix
    L <- unlist(eigen(J)$values) # eigenvalues of the Jacobian matrix
    maxreL <- max(Re(L)) # max real part of Jacobian matrix
    ct <- all(Neq>=0) & maxreL < (-0.005) & maxreL > - 0.1 # check for local stability and feasibility
    #ct <- maxreL < (-0.005) & maxreL > - 0.1

    if(as.numeric(ct)==1){
      k<-k+1
      community<-list(J,R_vec,A_mat,Neq,maxreL,as.matrix(new_spp_ls))
      names(community)<-c("jacobian_matrix","R","A_matrix","Neq","maxreL","tlvl")
      foodweb_out[[k]]<-community   
      setTxtProgressBar(pb, k)
      
    }else{
      k<-k
    } 
  } 
  close(pb)
  return(foodweb_out)
}

# assign_A_matrix_troph_comp<-function(motif.id,adjacent_matrix,new_spp_ls){
#   #A_matrix <- matrix(NA, 4, 4)
#   A_matrix <- matrix(NA,nrow = nrow(adjacent_matrix),ncol = ncol(adjacent_matrix))
#   
#   # assign the diagonal entries - intraspecific competiton
#   diagonal.elements <- diag(adjacent_matrix)
#   diagonal.elements[diagonal.elements==0] <- -0.1
#   diag(A_matrix) <- diagonal.elements
#   
#   # assign the interspecific competition between basal species
#   # how many basal species
#   basal.species.n <- length(which(diagonal.elements == -1)) 
#   # all entries are sampled from the uniform distribution [-0.5,0]
#   A_matrix_sub <- matrix(runif(basal.species.n^2,min=-0.5,max=0),basal.species.n,basal.species.n) 
#   # change the diagonal entries of basal species to -1
#   diag(A_matrix_sub) <--1  
#   A_matrix[1:basal.species.n,1:basal.species.n] <- A_matrix_sub
#   
#   # assign the consumer-resource interaction coefficient
#   for(i in (basal.species.n+1):nrow(adjacent_matrix)){ # basal species should not be a consumer
#     # how many consumer links
#     possible_links_from_species_i_to_j <- adjacent_matrix[1:(i-1), i]
#     N_possible_links_from_species_i_to_j <- length(which(possible_links_from_species_i_to_j==-1))
#     if(N_possible_links_from_species_i_to_j == 1){
#       possible_links_from_species_i_to_j[possible_links_from_species_i_to_j==-1] <- -0.5
#     }else{
#       consumer_interaction_value <- c(-0.4, rep(-0.1/(N_possible_links_from_species_i_to_j-1), N_possible_links_from_species_i_to_j-1))
#       #consumer_interaction_value <- c(-0.4, rep(-0.2/(N_possible_links_from_species_i_to_j-1), N_possible_links_from_species_i_to_j-1))
#       #consumer_interaction_value <- c(-0.4, runif(N_possible_links_from_species_i_to_j-1,-0.4,0))
#       permute_value <- sample(consumer_interaction_value, length(consumer_interaction_value))
#       possible_links_from_species_i_to_j[possible_links_from_species_i_to_j==-1] <- permute_value
#     }
#     A_matrix[1:(i-1), i] <- possible_links_from_species_i_to_j
#   }
#   
#   # assign the resource-consumer interaction coefficiens
#   # generate the energy conversion rate matrix first
#   conversion_matrix <- A_matrix
#   conversion_matrix[!is.na(conversion_matrix)] <- 1
#   conversion_matrix[is.na(conversion_matrix)] <- -0.2
#   
#   if(motif.id %in% c(7,10)){ # add lower influence in circumstances where omnivory
#     #conversion_matrix[4,1] <- -0.02
#     conversion_matrix[new_spp_ls[new_spp_ls$tlvl == 4,2],new_spp_ls[new_spp_ls$tlvl == 1,2]] <- -0.02
#     
#   }else if(motif.id %in% c(8,9)){
#     #conversion_matrix[3,1] <- -0.02
#     conversion_matrix[new_spp_ls[new_spp_ls$tlvl == 3,2],new_spp_ls[new_spp_ls$tlvl == 1,2]] <- -0.02
#     
#   }else if(motif.id %in% c(11,13)){
#     #conversion_matrix[3,1] <- -0.02
#     conversion_matrix[new_spp_ls[new_spp_ls$tlvl == 3,2],new_spp_ls[new_spp_ls$tlvl == 1,2]] <- -0.02
#     #conversion_matrix[4,1] <- -0.02
#     conversion_matrix[new_spp_ls[new_spp_ls$tlvl == 4,2],new_spp_ls[new_spp_ls$tlvl == 1,2]] <- -0.02
#     
#   }else if(motif.id == 12){
#     #conversion_matrix[4,1] <- -0.02
#     conversion_matrix[new_spp_ls[new_spp_ls$tlvl == 4,2],new_spp_ls[new_spp_ls$tlvl == 1,2]] <- -0.02
#     #conversion_matrix[4,2] <- -0.02
#     conversion_matrix[new_spp_ls[new_spp_ls$tlvl == 4,2],new_spp_ls[new_spp_ls$tlvl == 2,2]] <- -0.02
#     
#   }else if(motif.id == 6){
#     #conversion_matrix[4,2] <- -0.02
#     conversion_matrix[new_spp_ls[new_spp_ls$tlvl == 4,2],new_spp_ls[new_spp_ls$tlvl == 2,2]] <- -0.02
#     
#   }else{
#     conversion_matrix <- conversion_matrix
#   }
#   
#   # to obtain resource-consumer interaction
#   for(m in 1:ncol(adjacent_matrix)){
#     for(n in 1:ncol(adjacent_matrix)){
#       if(new_spp_ls[new_spp_ls$species == m,1] == new_spp_ls[new_spp_ls$species == n,1] & (is.na(A_matrix[m,n]) | A_matrix[m,n] == 0)){
#         A_matrix[m,n] <- runif(1,min=-0.5,max=0)
#       }else{
#         A_matrix[m,n] <- ifelse(is.na(A_matrix[m,n]),
#                                 conversion_matrix[m,n]*A_matrix[n,m],A_matrix[m,n])
#       }
#     }
#   }
#   return(A_matrix)
# }

#kk <- seqtime::generateA(N = 5, type = 'random',d=-0.5, min.strength = -0.5, max.strength = 0.5, c = 0.5,groups = c())
