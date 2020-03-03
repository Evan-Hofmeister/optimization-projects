NumAssets = 7

#
# Utility functions
#

# accumulation utility at final stage
myAccumulationUtilityFunction <- function(x,Target){if(x > 1.2 * Target) log(x) 
                                        else {if(x >= 0.8 * Target) log(x) - ((1.2*Target - x)^2)/1000000
                                              else log(x) - (Target^2)/6250000 }}
# 1-yr drawdown risk utility each year
myDrawdownUtilityFunction <- function(YearsBeforeRetire,nport){
        Drawdown[nport]*(20 - YearsBeforeRetire)*DrawdownAversion}


#
# Savings function for annual savings amounts
#
mySavingsFunction <- function(current_age, start_age = 22, start_salary = 60,
                              annual_salary_increase = 1.5, savings_rate = 0.14){
      if(current_age == RetirementAge) 0 
          else ((current_age - start_age) * annual_salary_increase + start_salary) * savings_rate}


#
# Look up to determine index in array corresponding
# to specific wealth level
#
# assumes wealth values in table are equally spaced
#
WealthIndex <- function(wealth,minwealth,maxwealth, wealthstep) 
                    {index = 1+round((wealth-minwealth)/wealthstep)
                          if(index < 1) 1 else 
                              if(index > 1+(maxwealth-minwealth)/wealthstep ) 1+round((maxwealth - minwealth)/wealthstep) 
                                else index}

##
## Inputs
##
##  You will probably want to change some of these
##

# allows me to store solution and avoid need
# to resolve as I figure out presentation
# materials and test different exhibits
# TRUE means solve, FALSE means read from files
Resolve = TRUE

z <- Sys.time()
MyFolder = "C:/Users/titan/Documents/CFRM 507"
WealthTarget = 126 * 0.55 / 0.025  ## 2,772,000
RetirementAge = 67
DrawdownAversion = 0.007

##
## Setup tree structure
##

# branching structure.  This tree only branches
# from the initial node
NumBranches = 1000
NumPeriods = 16
NumTrees <- 5

TreeStructure <- rep(1,NumPeriods)
TreeStructure[1] <- NumBranches

# Time horizon pattern.  Make sure it adds up to (retirement age - starting age)
# and that it is of same length as TreeStructure
TimeHorizonStructure <- rep(1,NumPeriods)

# wealth table range and interval
MinWealth = 100
MaxWealth = 6.0e3
WealthStep = 25

#
# return assumptions
#
AssetNames <- c("US Equity", "Int'l Equity", "Emer Mkt Equity", "Global RE", "Agg Fixed", "Hedge Fund","Cash")

ret_vec <- c(0.06, 0.061, 0.07, 0.052, 0.025, 0.052, 0.015)
stdev_vec <- c(0.191, 0.202, 0.268, 0.207, 0.068, 0.07, 0.058)
correl_mat <- rbind(c(1.00, 0.74, 0.67, 0.74, 0.13, 0.47, 0.02),
                    c(0.74, 1.00, 0.70, 0.78, 0.09, 0.46, 0.00),
                    c(0.67, 0.70, 1.00, 0.66, 0.07, 0.45, -0.03),
                    c(0.74, 0.78, 0.66, 1.00, 0.10, 0.37, -0.03),
                    c(0.13, 0.09, 0.07, 0.10, 1.00, 0.10, 0.10),
                    c(0.47, 0.46, 0.45, 0.37, 0.10, 1.00, 0.55),
                    c(0.02, 0.00, -0.03, -0.03, 0.10, 0.55, 1.00))
cov_mat <- diag(stdev_vec) %*% correl_mat %*% diag(stdev_vec)


print('Creating tree structure')

#
# Construct tree structure
#
FirstNodeAtStage = array(1,dim=c(length(TreeStructure)+1))
LastNodeAtStage = array(1,dim=c(length(TreeStructure)+1))
NodesAtStage = array(1,dim=c(length(TreeStructure)+1))
                                                    
TotalNodes = 1
                     
for(i in 2:(length(TreeStructure)+1))                 
{
  FirstNodeAtStage[i] = LastNodeAtStage[i-1] + 1
  NodesAtStage[i] = NodesAtStage[i-1] * TreeStructure[i-1]
  LastNodeAtStage[i] = FirstNodeAtStage[i] + NodesAtStage[i] - 1
  TotalNodes = TotalNodes + NodesAtStage[i]
}
ParentNode = array(0,dim=c(length(TotalNodes)))
FirstChild = array(2,dim=c(length(TotalNodes)))
LastChild = array(TreeStructure[1]+1,dim=c(length(TotalNodes)))
NodeStage = array(1,dim=c(length(TotalNodes)))
    
    
for(i in 1:LastNodeAtStage[length(TreeStructure)])
{
  for(j in FirstChild[i]:LastChild[i])
  {
    ParentNode[j] = i
    NodeStage[j] = NodeStage[ParentNode[j]] + 1
    if(NodeStage[j] <= length(TreeStructure))
    {
      FirstChild[j] = LastChild[j-1]+1
      LastChild[j] = FirstChild[j] + TreeStructure[NodeStage[j]] - 1
    }
  }
}


print('Creating scenarios')

#    
#  Construct Scenarios of Asset Returns                                     
#
library("MASS")

# reserve space to hold scenario values
AssetReturns = array(0,dim=c(NumTrees,TotalNodes,length(AssetNames)))

# draw samples from multi-variate lognormal distribution
for(ntree in 1:NumTrees)
{
  # set seed so you can replicate results, but choose a different seed for each tree
  set.seed(1234 * ntree)

  for(nperiod in 1:length(TreeStructure))
  {
    # determine statistics of associated normal distribution
    #  the normal variates correspond to the continuously compounded rates of return
    period_cc_stdev_vec = sqrt(log(1 + stdev_vec^2 / (1 + ret_vec)^2))
    period_cc_ret_vec = log(1+ret_vec) - period_cc_stdev_vec^2/2
    period_cc_stdev_vec = period_cc_stdev_vec * sqrt(TimeHorizonStructure[nperiod])
    period_cc_ret_vec = period_cc_ret_vec * TimeHorizonStructure[nperiod]
    period_cc_cov_mat = diag(period_cc_stdev_vec) %*% correl_mat %*% diag(period_cc_stdev_vec)
  
    # draw a sample from multi-variate normal
    AssetReturns[ntree, FirstNodeAtStage[nperiod+1]:LastNodeAtStage[nperiod+1],1:length(AssetNames)] = 
       mvrnorm(LastNodeAtStage[nperiod+1] - FirstNodeAtStage[nperiod+1] + 1, period_cc_ret_vec,period_cc_cov_mat)
  
    # convert to multi-variate lognormal
    AssetReturns[ntree,FirstNodeAtStage[nperiod+1]:LastNodeAtStage[nperiod+1],1:length(AssetNames)] = 
      exp(AssetReturns[ntree,FirstNodeAtStage[nperiod+1]:LastNodeAtStage[nperiod+1],1:length(AssetNames)]) - 1
    
    # diagnostics
    period_cc_ret_vec
    period_cc_stdev_vec
    colMeans(AssetReturns[ntree,FirstNodeAtStage[nperiod+1]:LastNodeAtStage[nperiod+1],])
  }
}



##
## Specify available portfolio mixes
##
NumPortfolios = 8
MyPortfolios = array(0,dim=c(NumPortfolios,length(AssetNames)))

MyPortfolios[1,1:dim(MyPortfolios)[2]] = c(0.30, 0.30, 0.10, 0.10, 0.10, 0.10, 0.00)
MyPortfolios[2,1:dim(MyPortfolios)[2]] = c(0.25, 0.25, 0.10, 0.10, 0.20, 0.10, 0.00)
MyPortfolios[3,1:dim(MyPortfolios)[2]] = c(0.25, 0.21, 0.08, 0.08, 0.30, 0.08, 0.00)
MyPortfolios[4,1:dim(MyPortfolios)[2]] = c(0.23, 0.19, 0.06, 0.06, 0.40, 0.06, 0.00)
MyPortfolios[5,1:dim(MyPortfolios)[2]] = c(0.21, 0.15, 0.04, 0.05, 0.45, 0.05, 0.05)
MyPortfolios[6,1:dim(MyPortfolios)[2]] = c(0.16, 0.12, 0.04, 0.04, 0.50, 0.04, 0.10)
MyPortfolios[7,1:dim(MyPortfolios)[2]] = c(0.13, 0.09, 0.02, 0.03, 0.55, 0.03, 0.15)
MyPortfolios[8,1:dim(MyPortfolios)[2]] = c(0.08, 0.06, 0.01, 0.02, 0.51, 0.02, 0.30)

#fee_vec <- c(0.0020, 0.0018, 0.0016, 0.0014, 0.0012, 0.0010, 0.0008, 0.0006)
#feeMat <- matrix(1,nrow=TotalNodes,ncol = 1) %*% fee_vec

# calculate node-by-node returns for portfolios
PortfolioReturns = array(0,dim=c(NumTrees,TotalNodes,NumPortfolios))
for(ntree in 1:NumTrees){
  AR = matrix(AssetReturns[ntree,1:TotalNodes,1:length(AssetNames)],nrow = TotalNodes, ncol = NumAssets)
  PortfolioReturns[ntree,1:TotalNodes,1:NumPortfolios] = AR %*% t(MyPortfolios) #- feeMat
}

## calculate drawdown risk
Drawdown <- apply(PortfolioReturns, 3, quantile, probs = 0.05)

print('Building wealth table')
print('Stage  Total  Time')
#
#  Backwards recursion
#

TotalWealths = 1 + (MaxWealth - MinWealth) / WealthStep
# arrays to hold key pieces of solution: 
#    stage (1 = age 67, 2 = age 66, ...)
#    wealth index
#    allocation mix of parent node
#    a separate array for 
#         1 = wealth level
#         2 = best decision (portfolio to hold)  no value for age 67
#         3 = expected utility at retirement
#SolutionArray = array(0,dim=c(length(TreeStructure)+1,TotalWealths,3))
# dimensions: # stages, # wealth levels, # portfolios
WealthArray = array(0,dim = c(length(TreeStructure)+1,TotalWealths))
DecisionArray = array(0,dim = c(length(TreeStructure)+1,TotalWealths))
EUtilityArray = array(0,dim = c(length(TreeStructure)+1,TotalWealths))

if(Resolve == TRUE)
{
  
  # Evaluate retirement age utility across wealth values
  # to initialize WealthArray and EUtilityArray
  #
  # This is the boundary condition and stage that is easy to solve
  for(nwealth in 1:TotalWealths)
{
  WealthArray[1,nwealth] = MinWealth + (nwealth-1) * WealthStep
  EUtilityArray[1,nwealth] = myAccumulationUtilityFunction(WealthArray[1,nwealth],WealthTarget)
}

# Run the recursion across each period using tree 1
for(nperiod in length(TreeStructure):1)
{
  for(nwealth in 1:TotalWealths)
  {
      WealthArray[length(TreeStructure) - nperiod + 2,nwealth] = MinWealth + (nwealth-1) * WealthStep
      MaxUtil = -1.0e20

      # check value for each mix and pick the best one
      for(nport2 in 1:NumPortfolios)
      {
        EUtil = 0
        for(nnode in FirstNodeAtStage[nperiod+1]:LastNodeAtStage[nperiod+1])
        {
          # wealth index at child node after returns and cash flow from parent node
          nwealthtplus1 = WealthIndex((WealthArray[length(TreeStructure) - nperiod + 2,nwealth]  + mySavingsFunction(RetirementAge + 1 - nperiod)) * (1+PortfolioReturns[1,nnode,nport2]),
                                      MinWealth,MaxWealth, WealthStep)

          EUtil = EUtil + EUtilityArray[length(TreeStructure) - nperiod+1,nwealthtplus1]
        }
        EUtil = EUtil / (LastNodeAtStage[nperiod+1] - FirstNodeAtStage[nperiod+1] + 1)
        EUtil = EUtil + myDrawdownUtilityFunction(length(TreeStructure) - nperiod + 2,nport2)
        if(EUtil > MaxUtil)
        {
          MaxUtil = EUtil
          BestChoice = nport2
        }
      
      DecisionArray[length(TreeStructure) - nperiod + 2,nwealth] = BestChoice
      EUtilityArray[length(TreeStructure) - nperiod + 2,nwealth] = MaxUtil
    }
  }
  print(c(nperiod, length(TreeStructure), Sys.time() - z))

}
} # end of if(Resolve == TRUE)

# read solution from files if it has already been solved
if(Resolve == FALSE)
{
  WealthArray <- read.table( paste(MyFolder,"SA_Wealth_2000.txt"), header=TRUE, sep="\t") 
  DecisionArray <- read.table(paste(MyFolder,"SA_Decision_2000.txt"), header=TRUE, sep="\t") 
  EUtilityArray <- read.table( paste(MyFolder,"SA_Utility_2000.txt"), header=TRUE, sep="\t")   
}

print('Creating client outcomes, tree 1')
print('Client       Mix    E[Wealth]  P{success} E[Util]   Actual Util')

##
##  Simulate results for each client on each tree
##

#  client characteristics
NumClients = 6
ClientAge = c(52, 55, 55, 57, 62, 65)
ClientBalance = c(500, 400, 900, 500, 1100, 950)
ClientBalance = c(550, 600, 950, 500, 1500, 1600)
ClientBalance = c(900, 700, 1000, 500, 1500, 1600)
ClientBalance = c(1000,2000, 500, 1500, 600, 2600 )

# values of initial decision and expected utility
ClientPortfolio = array(0, dim = c(6))
ClientEUtil = array(0, dim = c(6))

# reserve space to hold node-specific future values
ProjectWealth = array(0,dim=c(length(TotalNodes)))
ProjectDecision = array(0,dim=c(length(TotalNodes)))

##
## Loop over all nodes in each tree for each client to 
## evaluate client situation along scenario paths
##
for(nclient in 1:NumClients)
{
  for(ntree in 1:NumTrees)
    {
    windex = WealthIndex(ClientBalance[nclient]+mySavingsFunction(ClientAge[nclient]),MinWealth,MaxWealth, WealthStep)
    ProjectDecision[1] = DecisionArray[RetirementAge + 1 - ClientAge[nclient],windex]
    #ClientEUtil[nclient] = EUtilityArray[RetirementAge + 1 - ClientAge[nclient], windex]
    ProjectWealth[1] = ClientBalance[nclient] + mySavingsFunction(ClientAge[nclient])
    EUtil = 0
    for(nnode in 2:LastNodeAtStage[RetirementAge - ClientAge[nclient]])
    {
      windex = WealthIndex(ProjectWealth[ParentNode[nnode]], MinWealth, MaxWealth, WealthStep)
      ProjectWealth[nnode] = (ProjectWealth[ParentNode[nnode]]  + mySavingsFunction(ClientAge[nclient] + NodeStage[nnode] - 1)) *
                (1+PortfolioReturns[ntree,nnode,DecisionArray[RetirementAge + 3 - ClientAge[nclient] - NodeStage[nnode],windex]]) 
      windex = WealthIndex(ProjectWealth[nnode], MinWealth, MaxWealth, WealthStep)
      ProjectDecision[nnode] = DecisionArray[RetirementAge + 2 - ClientAge[nclient] - NodeStage[nnode],windex]
      EUtil <- EUtil + myDrawdownUtilityFunction(RetirementAge + 2 -ClientAge[nclient] - NodeStage[nnode],ProjectDecision[nnode])
    }
    AvgWealth = 0
    ProbGreaterThanTarget = 0
    for(nnode in FirstNodeAtStage[RetirementAge + 1 - ClientAge[nclient]]:LastNodeAtStage[RetirementAge + 1-ClientAge[nclient]])
    {
      windex = WealthIndex(ProjectWealth[ParentNode[nnode]], MinWealth, MaxWealth, WealthStep)
      ProjectWealth[nnode] = (ProjectWealth[ParentNode[nnode]]  + mySavingsFunction(ClientAge[nclient] + NodeStage[nnode] - 1)) *
        (1+PortfolioReturns[ntree,nnode,DecisionArray[RetirementAge + 3 - ClientAge[nclient] - NodeStage[nnode],windex]]) 
      AvgWealth = AvgWealth + ProjectWealth[nnode]
      if(ProjectWealth[nnode] > WealthTarget) ProbGreaterThanTarget = ProbGreaterThanTarget + 1
      EUtil = EUtil + myAccumulationUtilityFunction(ProjectWealth[nnode], WealthTarget)
    }
    AvgWealth = AvgWealth / (LastNodeAtStage[RetirementAge + 1 - ClientAge[nclient]] - FirstNodeAtStage[RetirementAge + 1-ClientAge[nclient]] + 1)
    ProbGreaterThanTarget = ProbGreaterThanTarget / (LastNodeAtStage[RetirementAge + 1 - ClientAge[nclient]] - FirstNodeAtStage[RetirementAge + 1-ClientAge[nclient]] + 1)
    EUtil = EUtil / (LastNodeAtStage[RetirementAge + 1 - ClientAge[nclient]] - FirstNodeAtStage[RetirementAge + 1-ClientAge[nclient]] + 1)
    EUtil = EUtil + myDrawdownUtilityFunction(RetirementAge + 1 - ClientAge[nclient],ProjectDecision[1])
    ClientEUtil[nclient] <- EUtil
    ClientPortfolio[nclient] <- ProjectDecision[1]
    # output to screen which client, recommended portfolio, average wealth, probability target exceeded, 
    #    expected utility based on solving, expected utility based on simulation
    print(c(nclient,ClientPortfolio[nclient],AvgWealth, ProbGreaterThanTarget, ClientEUtil[nclient],EUtil))
  }
}

#
#  Make some pretty pictures
#
myxlim <- c(500,3000)
#plot(WealthArray[2,],DecisionArray[2,],ylim = rev(c(1,NumPortfolios)),col = "blue", xlim = myxlim)
#plot(WealthArray[2,],DecisionArray[3,],ylim = rev(c(1,NumPortfolios)),col = "green", xlim = myxlim)
#plot(WealthArray[2,],DecisionArray[4,],ylim = rev(c(1,NumPortfolios)),col = "purple", xlim = myxlim)
#plot(WealthArray[2,],DecisionArray[6,],ylim = rev(c(1,NumPortfolios)),col = "red", xlim = myxlim)
#plot(WealthArray[2,],DecisionArray[8,],ylim = rev(c(1,NumPortfolios)),col = "yellow", xlim = myxlim)
#plot(WealthArray[2,],DecisionArray[10,],ylim = rev(c(1,NumPortfolios)),col = "orange", xlim = myxlim)
#plot(WealthArray[2,],DecisionArray[11,],ylim = rev(c(1,NumPortfolios)),col = "black", xlim = myxlim)

myTitle = paste("Drawdown aversion = ",DrawdownAversion,", NumBranches = ",NumBranches, sep="")
plot(WealthArray[2,],DecisionArray[2,], main = myTitle,
     ylim = rev(c(1,NumPortfolios)),col = "blue", xlim = myxlim)
lines(WealthArray[2,],DecisionArray[2,],col = "blue")
lines(WealthArray[2,],DecisionArray[3,],col = "green")
#lines(WealthArray[2,],DecisionArray[6,],col = "purple")
lines(WealthArray[2,],DecisionArray[6,],col = "red")
lines(WealthArray[2,],DecisionArray[11,],col = "yellow")
lines(WealthArray[2,],DecisionArray[13,],col = "orange")
lines(WealthArray[2,],DecisionArray[16,],col = "black")

WealthArray[2,]
DecisionArray[2,]

myTitle = paste("Drawdown aversion = ",DrawdownAversion,", NumBranches = ",NumBranches, sep="")
plot(WealthArray[2,],DecisionArray[2,], main = myTitle,
     ylim = rev(c(1,NumPortfolios)),col = "blue", xlim = myxlim)
lines(WealthArray[2,],DecisionArray[1,],col = "blue")
lines(WealthArray[2,],DecisionArray[2,],col = "blue")
lines(WealthArray[2,],DecisionArray[3,],col = "green")
lines(WealthArray[2,],DecisionArray[4,],col = "purple")
lines(WealthArray[2,],DecisionArray[5,],col = "red")
lines(WealthArray[2,],DecisionArray[6,],col = "yellow")
lines(WealthArray[2,],DecisionArray[7,],col = "orange")

myTitle = paste("Drawdown aversion = ",DrawdownAversion,", NumBranches = ",NumBranches, sep="")
plot(WealthArray[2,],DecisionArray[2,], main = myTitle,
     ylim = rev(c(1,NumPortfolios)),col = "blue", xlim = myxlim)
lines(WealthArray[4,],DecisionArray[1,],col = "blue")
lines(WealthArray[4,],DecisionArray[2,],col = "blue")
lines(WealthArray[4,],DecisionArray[3,],col = "green")
lines(WealthArray[4,],DecisionArray[4,],col = "purple")
lines(WealthArray[4,],DecisionArray[5,],col = "red")
lines(WealthArray[4,],DecisionArray[6,],col = "yellow")
lines(WealthArray[4,],DecisionArray[7,],col = "orange")




# store solution so it can be read back in later
if(Resolve == TRUE)
{
  write.table(WealthArray[,], paste(MyFolder,"SA_Wealth_2000.txt"), sep="\t") 
  write.table(DecisionArray[,], paste(MyFolder,"SA_Decision_2000.txt"), sep="\t") 
  write.table(EUtilityArray[,], paste(MyFolder,"SA_Utility_2000.txt"), sep="\t") 
}
SolutionTime = Sys.time() - z
print(SolutionTime)