library(ape)


mytree <- read.tree("/path/to/rooted/tree") 

node <- c(
        getMRCA(mytree, tip = c("Oryziaslatipes","Takifugurubripes") ),
        getMRCA(mytree, tip = c("Latescalcarifer","Takifugurubripes") ),
        getMRCA(mytree, tip = c("Larimichthyscrocea","Takifugurubripes") ),
        getMRCA(mytree, tip = c("Oreochromisniloticus","Takifugurubripes") )
        )

# Ages taken from TIMETREE (http://timetree.org) for each corresponding pair of species
age.min <- c(
            104,
            94,
            104,
            104
        )

age.max <- c(
            145,
            115,
            145,
            145
        )

soft.bounds <- c(
            FALSE,
            FALSE,
            FALSE,
            FALSE
        )

mycalibration <- data.frame(node, age.min, age.max, soft.bounds)

mytimetree <- chronos(mytree, lambda = 1, model = "discrete", calibration = mycalibration, control = chronos.control() )

write.tree(mytimetree,file="/path/to/save/the/CALIBRATED.tree")

