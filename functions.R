
# Ising Resting State Distibution -----------------------------------------------

ising.dist <- function(omega, mu, beta)
{
  n = length(mu) # number of magnets
  x = rep(0,n) # set binary class of magnets (0/1)

  foo <- function(x)
  {
    if(x==-1){y = c(x,0)}
    if(x==0){y = c(x,1)}
    if(x==1){y = c(-1,x)}
    return(list(y))
  }

  x = expand.grid(sapply(x,foo)) # create all possible patterns

  boltzmann <- function(x,b=1){exp(-b*x)} # boltzmann function
  hamiltonian <- function(x,mu,omega){(1/2)*(-t(x)%*%omega%*%x) - t(x)%*%mu} # hamiltonian function

  E = apply(x,1,hamiltonian,mu=mu,omega=omega) # calculate energy for each configuration
  B = boltzmann(E,b = beta) # apply boltzmann function to energies
  P = B/sum(B) # calculate stationary probabilities
  y = cbind(P,x) # combine probabilities with patterns
  colnames(y)[1] = c("Probability")
  return(y)
}


# Expected Choice Behaviour -----------------------------------------------

ecb <- function(omega, mu, beta, vareps, plot = FALSE, n.cue = 1)
{

  ## INPUT
  # omega = the pairwise interaction matrix
  # mu = the base appeal for the cue and alternatives
  # beta = inverse temperature
  # vareps = lateral inhinbition

  # load dependencies
  require("plyr")

  n.cue = n.cue # number of cues
  n.node = ncol(omega) # detect number of nodes
  n.alt = n.node-n.cue # set number of alternatives

  # create a plot of the choice setup
  if(plot)
  {
    dev.off()
    require("qgraph")
    cue.lay <- seq(0, (n.cue-1)*2,by=2)
    cue.lay <- cue.lay-mean(cue.lay)

    alt.lay <- seq(0, (n.alt-1)*2,by=2)
    alt.lay <- alt.lay-mean(alt.lay)

    adj.lay <- matrix(c(cue.lay,alt.lay,rep(2,n.cue),rep(0,n.alt)),n.cue + n.alt,2)

    curve.setup = qgraph(omega,layout=adj.lay,DoNotPlot=T)
    graph.curve = apply(cbind(unlist(curve.setup$Edgelist[1]),unlist(curve.setup$Edgelist[2])),1,function(x)all(x>n.cue))
    plot.maximum =  max(abs(c(as.vector(omega),as.vector(mu))))*2
    qgraph(omega,layout=adj.lay,curve=-as.numeric(graph.curve)*3,mar=rep(10,4),labels=paste(colnames(omega),"\n mu = ",round(mu,2)),maximum = plot.maximum,weighted=T,negDashed=T,palette="gray")
  }

  # resting state distribution full structure
  full.x = ising.dist(omega,mu,beta)

  # resting state distribution alternatives only
  alt.columns = names(full.x)[-c(1:(n.cue+1))]
  alt.x = plyr::ddply(full.x, alt.columns, plyr::summarize, prob=sum(Probability))
  alt.start.prob = alt.x[,ncol(alt.x)]

  # all possible alternative patterns
  alt.x <- as.matrix(alt.x[,-ncol(alt.x)])
  colnames(alt.x) = colnames(omega)[(n.cue+1):n.node]

  # add lateral inhibition to weight.matrix
  omega.vareps <- omega
  omega.vareps[-c(1:n.cue),-c(1:n.cue)] = omega.vareps[-c(1:n.cue),-c(1:n.cue)]-vareps
  diag(omega.vareps) = 0

  # metropolis function - single-spin-flip dynamics (glauber dynamics)
  ap.gb <- function(x,alt.x,n.cue,omega,mu,beta,n.alt,no.iter)
  {

    hamiltonian <- function(x, omega, mu, beta){-beta*(((1/2)*(-t(x)%*%omega%*%x)) - t(x)%*%mu)}

    ## input
    i = x[1] # index to which alternative configuration is the current state
    j = x[2] # index to which alternative is selected for flipping
    k = x[3] # information for progressbar

    setTxtProgressBar(pb,k)

    x.c = alt.x[i,] # grab current configuration based on i
    x.p = x.c # make new variable for proposal configuration

    x.p[j] = (x.p[j]-1)^2 # flip value of selected variable

    # if the proposal is accepted which alternative configuration would we get to.
    to = which(apply(alt.x,1,function(y) all.equal(y,x.p)=="TRUE"))

    xf.c = c(rep(1,n.cue),x.c) # add cue states (all = 1) to current state
    xf.p = c(rep(1,n.cue),x.p) # add cue states (all = 1) to proposal state

    # acceptance probability
    ap <- min(1,exp(hamiltonian(xf.p, omega, mu,beta)-hamiltonian(xf.c, omega, mu,beta)))

    # correct the acceptance probability for the selection probability.
    ap <-  ap/n.alt

    # return alternative configuration if proposal is accepted, and acceptance probability of proposal.
    return(c(to,ap))
  }

  # apply function acceptance.probability to all combinations of configurations, alternatives, and proposals.
  iter <- expand.grid(1:nrow(alt.x),1:n.alt)
  iter <- cbind(iter,1:nrow(iter))
  no.iter = nrow(iter)

  cat("Calculating Transition Probabilities \n")
  pb = txtProgressBar(min = 0, max = no.iter, initial = 0, style = 3)
  out = cbind(iter[,1],t(apply(iter,1,ap.gb,alt.x,n.cue,omega.vareps,mu,beta,n.alt,no.iter)))

  # create the transition matrix with acceptance probabilities,
  trans.mat <- matrix(0,nrow(alt.x),nrow(alt.x))
  trans.mat[as.matrix(out[,1:2])] = out[,3]

  cat("\n")

  colnames(trans.mat) = rownames(trans.mat) = apply(alt.x,1,function(x)paste0(x,collapse="")) # add acceptance probabilities to the transition matrix
  diag(trans.mat) = diag(trans.mat) + (1-rowSums(trans.mat)) # incorporate rejection probabilities into the transition matrix diagonal

  absorbing.states <- which(apply(alt.x,1,function(x)length(which(x==1)))==1) # in which configurations is only one alternative active
  transitive.states <- which(apply(alt.x,1,function(x)length(which(x==1)))!=1) # in which configurations are none or more than one alternative active

  ## create transition matrix for first-stop paradigm, by setting transition probabities to 1 on the diagonal for absorbing states.
  first.stop.mat <- trans.mat
  first.stop.mat[absorbing.states,] = 0
  diag(first.stop.mat) = diag(first.stop.mat)+(1-rowSums(first.stop.mat))
  diag(first.stop.mat)[absorbing.states] = 1

  ## Extract canonical parts of transition matrix

  # matrix with transition probabilities from transitive to transitive states
  P = first.stop.mat[transitive.states,transitive.states]
  # matrix with transition probabilities from transitive to absorbing states
  R = first.stop.mat[transitive.states,absorbing.states]
  # matrix with transition probabilities from absorbing to absorbing states (i.e., identity matrix)
  I = first.stop.mat[absorbing.states,absorbing.states]

  zt = alt.start.prob[transitive.states] # starting probabilities for transitive configurations
  za = alt.start.prob[absorbing.states] # starting probabilities for absorbing configurations

  # expected response proportions for each of the alternatives
  ER = zt %*% Matrix::solve(diag(nrow(P))-P) %*% R
  ER = alt.start.prob[absorbing.states]+ER # correct for probabilities of the process starting in an absorbing state
  colnames(ER) = paste0("ALT",apply(alt.x[absorbing.states,],1,function(x)which(x==1)))
  ER = ER[,sort.int(colnames(ER),index.return = T)$ix] # sort the output by colname

  # expected response time (number of iterations) for process is absorbed at each of the alternatives
  ET = zt %*% Matrix::solve(diag(nrow(P))-P) %*% Matrix::solve(diag(nrow(P))-P) %*% R
  colnames(ET) = paste0("ALT",apply(alt.x[absorbing.states,],1,function(x)which(x==1)))
  ET = ET[,sort.int(colnames(ET),index.return = T)$ix]
  ET = ET/ER

  # calculate conditional distribution of chain as a function of the base appeal and connectivity matrix
  #ccp <- exp(beta*(mu[-1] + omega[-c(1:n.cue),-c(1:n.cue)]))/sum(exp(beta*(mu[-1] + omega[-1,1])))

  if(n.cue>1)
  {
    ccp <- exp(beta*(mu[-c(1:n.cue)] + rowSums(omega[-c(1:n.cue),c(1:n.cue)])))/sum(exp(beta*(mu[-c(1:n.cue)] + rowSums(omega[-c(1:n.cue),c(1:n.cue)]))))
  } else {
    ccp <- exp(beta*(mu[-1] + omega[-1,1]))/sum(exp(beta*(mu[-1] + omega[-1,1])))
  }

  # combine expected responses and response times for 1st top condition, and response probabilities for convergence of chain as a function of parameters and eigenvalue decomposition of the transition matrix.
  stat.dist <- rbind(ccp,ER,ET)
  colnames(stat.dist) = colnames(omega)[(n.cue+1):n.node]
  rownames(stat.dist) = c("conditional choice probabilities","first-stop choice probabilities","first-stop expected time")

  output <- list(stat.dist)
  names(output) <- c("results")

  print(round(stat.dist,2))
  return(output)
}


