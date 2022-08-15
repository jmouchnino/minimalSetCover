# Projet Problème de couverture mnimal par ensembles

# Soit U l'ensemble constitué des n premiers entiers {1, 2, ..., n}.
# S un ensemble de parties de U.
# On cherche un sous-emsemble T de S tel que 
#l'union des éléments de T est égal à U
library(magrittr)

sample_subsets <- function(n_sets, n) {
  # Prenons par exemple  S un ensemble de n parties de U, 
  # chaque partie de taille racine carré de n
  nbParties = n_sets
  lambda = sqrt(n)
  CardPartie = rpois(nbParties, lambda) 
  
  # Construction de S
  S <- lapply(seq_len(nbParties), \(i) sample(1:n, CardPartie[i]))
  return(S)
}

#### Première méthode ########### 
min_set_cover_greed <- function(n, nmax = 20000, beta0 = 0.05, d) {
  U = 1:n
  # Prenons par exemple  S un ensemble de n parties de U, 
  # chaque partie de taille racine carré de n
  nbParties = n 
  lambda = sqrt(n)
  CardPartie = rpois(nbParties, lambda) 
  # hist(CardPartie, freq = FALSE)
  # Construction de S
  S = vector("list", length = nbParties)
  for(i in 1:nbParties){
    nouvelPartie = sample(1:n, CardPartie[i] , replace = FALSE)
    S[[i]] = nouvelPartie
  }
  
  # Vérifions tout d'abord qu'on peut recouvrir tous les
  # éléments de U avec l'ensemble S
  
  E = unique(unlist(S)) # éléments couverts par G
  nonCouvertS = length(setdiff(U,E))
  nonCouvertS
  
  # S couvre bien tous les éléments de U
  
  # Fonction à minimiser
  h= function(G){
    E = unique(unlist(G)) # éléments couverts par G
    return(length(setdiff(U,E))) # setdiff renvoie les éléments qui ne sont pas dans E
  }
  
  # essayons de déterminer cardTmax
  cardTmax = length(S) # initialisation
  cardTmax
  cardTmin = 1
  difference = 10
  while(cardTmax - cardTmin > 1){
    cardT = ((cardTmin + cardTmax)/2 )%>% floor
    print(cardT)
    G=sample(S, size = cardT, replace=F)
    transition = function(Gin){
      Gout = Gin
      coord_choisie = sample (1:cardT, size = 1)
      Gout[coord_choisie]= sample(S, size = 1)
      return(c(Gout, coord_choisie))
    }
    # stockage de H
    H = rep(NA, nmax)
    H [1]=h(G)
    H_best = H[1]
    G_best = G
    for(k in 2:nmax){
      beta =beta0*sqrt(k)
      toto =transition(G)
      # calcul de delta
      Gprime = toto[1:cardT]
      K0 = toto[cardT+1]
      delta = h(Gprime)- h(G)
      # test
      u =runif(1)
      if ( u < exp(-beta*delta)) { #faut il mettre un moins devant beta ?
        G = Gprime
        H[[k]] =H[[k-1]]+ delta
      }
      else {
        H[[k]] = H[[k-1]]
      }
      if(H[[k]]< H_best){
        H_best =H[[k]]
        G_best = G
      }
      if(H_best ==0){ # sortir de la boucle si H a atteint 0
        Gpluspetit = G_best
        break
      } 
    }
    cat("H_best =", H_best, " ")
    cat("card_Gbest =", sum(G_best), " ")
    cardT = ((cardTmin + cardTmax)/2 )%>% floor
    if (H_best ==0) {
      cardTmax = cardT
    }
    if (H_best != 0){
      cardTmin = cardT
    }
  }
  cardTmax
  cardTmin
  cardT
  # plot(1:nmax, H, type = "line")
  Gpluspetit
  
}
#### Deuxième méthode ########### 
min_set_cover_penal <- function(
    n, 
    n_sets = n,  
    n_iter = 10000, 
    gamma = 1.1, 
    beta0 = 1, 
    proportion = runif(1),
    G= sample(c(0,1), size = n_sets, prob = c(1-proportion, proportion), replace=T)) {
  
  U = 1:n
  nbParties = n_sets
  S <- sample_subsets(n_sets, n)
  
  # On considère que G est un vecteur de longeur nbParties (length(S)).
  # On met un 1 si le sous-ensemble est présent et
  # et un 0 sinon
  
  gamma = 1.9 # coefficient pour la forme lagrangienne
  
  # Prendre un nombre légèrement supérieur à 1
  
  # Fonction à minimiser
  h2 = function(G){
    nonCouv <- S[G==1] %>% unlist %>% unique %>%
      setdiff(U, .) %>% length() 
    return(sum(G) + gamma * nonCouv)
  }
  
  transition2 = function(G_in){
    element_choisie <- sample(1:nbParties, size = 1) 
    G_in[element_choisie] <- 1 - G_in[element_choisie] # si 1 changer en 0 et inversement
    return(G_in)
  }
  
  # Intitialisation
  # proportion = runif(1) # proportion de 1
  # G=sample(c(0,1), size = nbParties, prob = c(1-proportion, proportion),replace=T)
  
  # beta
  nmax = n_iter
  H = rep(NA, nmax)
  H[1]=h2(G)
  H_best = H[1]
  G_best = G
  for(k in 2:nmax){
    if(k %% 100 == 0) cat(k, "card_Gbest =", sum(G_best), " ")
    
    beta =beta0*sqrt(k)
    # calcul de delta
    Gprime = transition2(G)
    delta = h2(Gprime)- h2(G)
    # test 
    if (runif(1)< exp(-beta * delta)) {
      G = Gprime
      H[k] = H[k-1] + delta
    } else {
      H[k] = H[k-1]
    }
    if(H[k]< H_best){
      H_best = H[k]
      G_best = G
    }
  }
  # G_best
  cardinalGminimlal = sum(G_best)
  # cardinalGminimlal # le cardinal de l'ensemble mnimal qui recouvre U
  iteration = 15000
  Gensemble <- S[G==1]
  # Gminimal
  E = unique(unlist(Gensemble))# éléments couverts par G
  
  # On vérifie que le nombre d'entier non couverts est égal à 0
  noncouv = length(setdiff(U,E)) # nombre d'éléments non couverts
  if (noncouv == 0) cat("\nAll the elements are covered ") else cat("\nSome elements are not covered")
  cat("\nCardinal of the minimal set cover :", sum(G), "\n")
  return(list(set_parties = S, ensemble = Gensemble, binary = G, couvert = unique(unlist(Gensemble))))
}
min_set_cover_penal(n = 10000, n_sets = 500, beta0 = 0.3,  n_iter = 10000, proportion = 0.001) %>% invisible()
