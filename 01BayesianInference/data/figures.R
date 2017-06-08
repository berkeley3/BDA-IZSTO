
kc8084 <- read.delim('gd80to84.txt')
names(kc8084) <- names(kc8084)[-grep('x', tolower(names(kc8084)))]

kc8589 <- read.delim('gd85to89.txt')

kc <- data.frame(state=kc8589$state, 
                 county=kc8589$county, 
                 aakc=kc8084$aadc + kc8589$aadc
                 #aakc=10^5*(kc8084$dc + kc8589$dc)/(kc8084$pop + kc8589$pop)
)

q90 <- with(kc, quantile(aakc, probs = 0.9, na.rm=T))
kc90 <- subset(kc, aakc>q90)

counties.90 <- paste0(tolower(kc90$state),',',tolower(kc90$county))


library(maps)
map('state')
map('county', counties.90, col='black', fill=T, add=T, main='label')
title('Highest kidney cancer death rates')

q10 <- with(kc, quantile(aakc, probs = 0.1, na.rm=T))
kc10 <- subset(kc, aakc<q10)

counties.10 <- paste0(tolower(kc10$state),',',tolower(kc10$county))


library(maps)
map('state')
map('county', counties.10, col='black', fill=T, add=T)
title('Lowest kidney cancer death rates')
