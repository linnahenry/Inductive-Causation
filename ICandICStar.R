library(igraph)

realmat <- matrix(0,12,12)

realmat[1,2] <- 1
realmat[2,3] <- 1
realmat[3,4] <- 1
realmat[4,5] <- 1
realmat[5,6] <- 1


realmat[7,2] <- 1
realmat[8,3] <- 1
realmat[9,4] <- 1
realmat[10,5] <- 1
realmat[11,6] <- 1


realmat[1,8] <- 1
realmat[2,9] <- 1
realmat[3,10] <- 1
realmat[4,11] <- 1
realmat[5,12] <- 1


realmat[7,8] <- 1
realmat[8,9] <- 1
realmat[9,10] <- 1
realmat[10,11] <- 1
realmat[11,12] <- 1

nn = 1000

datamat <- matrix(0, nn, 12) #n = 100
errormat <- matrix(0, nn, 12) #n = 100 #p=0.05, 0.2
datamat[,1] <- rbinom(nn, 1, 0.25) #n=100 
datamat[,7] <- rbinom(nn, 1, 0.25) #n=100 

prob = 0.2
for(i in 1:12){
	errormat[,i] <- rbinom(nn, 1, prob)	
}	

x <- datamat[,1] | datamat[,7]
datamat[(x == T & errormat[,1] == 1),2] <- 0
datamat[(x == T & errormat[,1] == 0),2] <- 1
datamat[(x == F & errormat[,1] == 1),2] <- 1
datamat[(x == F & errormat[,1] == 0),2] <- 0

for(i in 2:6){
	j = i+6
	x <- datamat[,(i-1)] | datamat[,(j-1)]
	datamat[(x == T & errormat[,(i-1)] == 1),i] <- 0
	datamat[(x == T & errormat[,(i-1)] == 0),i] <- 1
	datamat[(x == F & errormat[,(i-1)] == 1),i] <- 1
	datamat[(x == F & errormat[,(i-1)] == 0),i] <- 0
	
	datamat[(x == T & errormat[,(j-1)] == 1),j] <- 0
	datamat[(x == T & errormat[,(j-1)] == 0),j] <- 1
	datamat[(x == F & errormat[,(j-1)] == 1),j] <- 1
	datamat[(x == F & errormat[,(j-1)] == 0),j] <- 0
}	
	

respmat <- matrix(0, 12, 12)
sepsets <- matrix(0, 12, 12)
sepsets <- as.data.frame(sepsets)
for(i in 1:12){
	for(j in i:12){
		#print(i)
		#print(j)
		tmat <- matrix(0, 2, 2)
		tmat[1,1] <- sum(datamat[,i] == 0 & datamat[,j] == 0)
		tmat[1,2] <- sum(datamat[,i] == 1 & datamat[,j] == 0)
		tmat[2,1] <- sum(datamat[,i] == 0 & datamat[,j] == 1)
		tmat[2,2] <- sum(datamat[,i] == 1 & datamat[,j] == 1)
		pv <- chisq.test(tmat)$p.val
		
		if(pv < 0.05){ #they are associated
			xv <- 1:12
			vtolook <- xv[-i]
			vtolook <- vtolook[vtolook!=j]
			res <- matrix(0, 10, 10)
			for(k in 1:9){
				for(m in k:10){
					stratvar <- datamat[,vtolook[k]]
					stratvar2 <- datamat[,vtolook[m]]
					#zeros <- cbind(datamat[stratvar == 0,i], datamat[stratvar==0, j])
					#firsts <- cbind(datamat[stratvar == 1,i], datamat[stratvar==1, j])
					zeros <- cbind(datamat[stratvar==0 & stratvar2==0,i], datamat[stratvar==0 & stratvar2==0,j])
					zeroone <- cbind(datamat[stratvar==0&stratvar2==1,i], datamat[stratvar==0 & stratvar2==1,j])
					onezero <- cbind(datamat[stratvar==1&stratvar2==0,i], datamat[stratvar==1 & stratvar2==0,j])
					firsts <- cbind(datamat[stratvar==1 & stratvar2==1,i], datamat[stratvar==1 & stratvar2==1,j])
				
					zmat <- matrix(0, 2, 2)
					zmat[1,1] <- sum(zeros[,1] == 0 & zeros[,2] == 0) + 1
					zmat[1,2] <- sum(zeros[,1] == 1 & zeros[,2] == 0) + 1
					zmat[2,1] <- sum(zeros[,1] == 0 & zeros[,2] == 1) + 1
					zmat[2,2] <- sum(zeros[,1] == 1 & zeros[,2] == 1) + 1
					zpv <- chisq.test(zmat)$p.val
				
					zomat <- matrix(0,2,2)
					zomat[1,1] <- sum(zeroone[,1] == 0 & zeroone[,2] == 0) + 1
					zomat[1,2] <- sum(zeroone[,1] == 1 & zeroone[,2] == 0) + 1
					zomat[2,1] <- sum(zeroone[,1] == 0 & zeroone[,2] == 1) + 1
					zomat[2,2] <- sum(zeroone[,1] == 1 & zeroone[,2] == 1) + 1
					zopv <- chisq.test(zomat)$p.val
					
					ozmat <- matrix(0,2,2)
					ozmat[1,1] <- sum(onezero[,1] == 0 & onezero[,2] == 0) + 1
					ozmat[1,2] <- sum(onezero[,1] == 1 & onezero[,2] == 0) + 1
					ozmat[2,1] <- sum(onezero[,1] == 0 & onezero[,2] == 1) + 1
					ozmat[2,2] <- sum(onezero[,1] == 1 & onezero[,2] == 1) + 1
					ozpv <- chisq.test(ozmat)$p.val
				
					fmat <- matrix(0,2,2)
					fmat[1,1] <- sum(firsts[,1] == 0 & firsts[,2] == 0) + 1
					fmat[1,2] <- sum(firsts[,1] == 1 & firsts[,2] == 0) + 1
					fmat[2,1] <- sum(firsts[,1] == 0 & firsts[,2] == 1) + 1
					fmat[2,2] <- sum(firsts[,1] == 1 & firsts[,2] == 1) + 1
					fpv <- chisq.test(fmat)$p.val
				
					pvs <- c(zpv, fpv, zopv, ozpv)
					res[k,m] <- sum(pvs <= 0.05)/4	

					if(res[k,m] == 0){ #then vtolook[k] and vtolook[m] are sep sets for i and j
						if(sepsets[i,j] == 0 | sepsets[i,j] == "0"){
							sepsets[i,j] <- paste(vtolook[k], ",", vtolook[m], "|", sep = "")
						}
						else{
							sepsets[i,j] <- paste(sepsets[i,j], vtolook[k], ",", vtolook[m], "|", sep = "")
						}
					}
					
					#if(zpv <= 0.05 & fpv <= 0.05){ #are associated
					#	res[k] <- 1
					#}	
					#else if(zpv > 0.05 & fpv > 0.05){ #are not associated
					#	res[k] <- 0	
					#}
					#else{
					#	res[k] <- 0.5	
					#}		
				}
				
			}	
			#respmat[i,j] <- sum(res)
			if(sum(res == 0) > 46){ #then zero at [i,j] else 1
				respmat[i,j] <- 0
			}else{
				respmat[i,j] <- 1
			}
	
				
		}else{ #not associated
			#do nothing
		}	
			
	}
}	

respmat
for(i in 1:12){
	if(respmat[i,i] == 1){
		respmat[i,i] <- 0	
	}		
}	
sepsets

utrealmat <- realmat
utrealmat[2,7] <- 1
utrealmat[3,8] <- 1
utrealmat[4,9] <- 1
utrealmat[5,10] <- 1
utrealmat[6,11] <- 1
utrealmat[7,2] <- 0
utrealmat[8,3] <- 0
utrealmat[9,4] <- 0
utrealmat[10,5] <- 0
utrealmat[11,6] <- 0

#h0: independent
#hA: not independent

#the higher the number in respmat, the more associated they are.


fullrespmat <- respmat
for(i in 1:11){
	for(j in (i+1):12){
		if(respmat[i,j] == 1){
			fullrespmat[j,i] <- 1
		}
	}		
}

#for each pair of non adjacent variables a and b with common neighbor c, check if c in Sab
#if it is, then continue
#if it is not, then point arrows a -> c <- b
#ufullrespmat <- fullrespmat
#index <- 1:12
#for(i in 1:11){
#	for(j in (i+1):12){
#		if(respmat[i,j] == 0){
#			ineighb <- index[fullrespmat[i,] ==1]
#			jneighb <- index[fullrespmat[j,] ==1]
#			#if(length(ineighb) <= length(jneighb){
#				if((sum(ineighb %in% jneighb) >= 2 | sum(jneighb %in% ineighb) >= 2) & (nchar(sepsets[i,j]) == 1)){
#					next
#				}
#				if(sum(ineighb %in% jneighb) >= 2){
#					ss <- ineighb[ineighb %in% jneighb]
#					rss <- sepsets[i,j]
#					smat <- matrix(0, 2, as.integer(nchar(rss)/4))
#					start <- 0
#					ki <- 1
#					for(k in 1:(as.integer(nchar(rss)))){
#						if(substr(rss, k, k) == ","){
#							smat[1,ki] <- as.numeric(substr(rss, (start+1), (k-1)))
#							start <- k
#						}
#						if(substr(rss, k, k) == "|"){
#							smat[2,ki] <- as.numeric(substr(rss, (start+1), (k-1)))
#							start <- k
#							ki <- ki+1
#						}
#			
#						#ki <- (k-1)*4 + 1
#						#smat[1,k] <- as.numeric(substr(rss, ki, ki))
#						#smat[2,k] <- as.numeric(substr(rss, (ki+2), (ki+2)))
#					}
#					for(k in 1:(length(ss)-1)){
#						for(m in (k+1):length(ss)){
#							if(sum(smat[1,] == ss[k] & smat[2,] == ss[m]) >= 1){
#								#do nothing
#							}else{
#								ufullrespmat[i, ss[k]] <- 1
#								ufullrespmat[i, ss[m]] <- 1
#								ufullrespmat[j, ss[k]] <- 1
#								ufullrespmat[j, ss[m]] <- 1
#								ufullrespmat[ss[k], i] <- 0
#								ufullrespmat[ss[m], i] <- 0
#								ufullrespmat[ss[k], j] <- 0
#								ufullrespmat[ss[m], j] <- 0
#							}
#						}
#					}		
#				}#otherwise, no common neighbors.
#			#}else{
#			#	if(sum(jneighb %in% ineighb) >= 2){
#			#	
#			#	}
#			#}
#		}
#	}
#}

#redo, check c's one at a time
ufullrespmat <- fullrespmat
index <- 1:12
for(i in 1:11){
	for(j in (i+1):12){
		if(respmat[i,j] == 0){
			ineighb <- index[fullrespmat[i,] ==1]
			jneighb <- index[fullrespmat[j,] ==1]
			ss <- ineighb[ineighb %in% jneighb]
			if(length(ss) == 0){
				next	
			}	
			rss <- sepsets[i,j]
			smat <- matrix(0, 2, as.integer(nchar(rss)/4))
			start <- 0
			ki <- 1
			for(k in 1:(as.integer(nchar(rss)))){
				if(substr(rss, k, k) == ","){
					smat[1,ki] <- as.numeric(substr(rss, (start+1), (k-1)))
					start <- k
				}
				if(substr(rss, k, k) == "|"){
					smat[2,ki] <- as.numeric(substr(rss, (start+1), (k-1)))
					start <- k
					ki <- ki+1
				}
	
				#ki <- (k-1)*4 + 1
				#smat[1,k] <- as.numeric(substr(rss, ki, ki))
				#smat[2,k] <- as.numeric(substr(rss, (ki+2), (ki+2)))
			}
			for(k in 1:length(ss)){
				if(!(ss[k] %in% smat)){
					ufullrespmat[i,ss[k]] <- 1
					ufullrespmat[j,ss[k]] <- 1
					ufullrespmat[ss[k],i] <- 0
					ufullrespmat[ss[k],j] <- 0		
				}	 	
			}	
			
		}
	}		
}

#r1
ufrm <- ufullrespmat
for(i in 1:11){
	for(j in (i+1):12){
		if(ufrm[i,j] == ufrm[j,i] & ufrm[i,j] == 1){
			print(paste(i,j))
			ipar <- index[ufrm[,i]==1]
			if(length(ipar) ==0){
				next
			}
			for(k in 1:length(ipar)){
				if(ipar[k] == j){
					next
				}
				if(ufrm[i, ipar[k]] == 0){ #directed edge from parent to i
					ufrm[i,j] <- 1
					ufrm[j,i] <- 0
				}
			}
			jpar <- index[ufrm[,j] ==1]
			if(length(jpar) ==0){
				next
			}
			for(k in 1:length(jpar)){
				if(jpar[k] == i){
					next
				}
				if(ufrm[j, jpar[k]] == 0){ #directed edge from parent to j
					ufrm[j,i] <- 1
					ufrm[i,j] <- 0
				}
			}
		}
	}
}
#I find that after R1, all edges are directed in every case.

#prints undirected edges
for(i in 1:11){
	for(j in (i+1):12){
		if(ufrm[i,j] == ufrm[j,i] & ufrm[i,j] == 1){
			print(paste(i,j))
		}
	}
}
#no edges printed (when nn=1000)!

g1 <- graph.adjacency(realmat, mode = "directed") #what the graph *should be*
#plot(g1)

rownames(ufrm) <- c("x0", "x1", "x2", "x3", "x4", "x5", "y0", "y1", "y2", "y3", "y4", "y5")
colnames(ufrm) <- c("x0", "x1", "x2", "x3", "x4", "x5", "y0", "y1", "y2", "y3", "y4", "y5")

g2 <- graph.adjacency(ufrm, mode = "directed")
plot(g2, main = paste("graph of n=", nn, "and p=", prob)) #end result after IC algo







#IC* algo
ufrm2 <- fullrespmat
index <- 1:12
for(i in 1:11){
	for(j in (i+1):12){
		if(respmat[i,j] == 0){
			ineighb <- index[fullrespmat[i,] ==1]
			jneighb <- index[fullrespmat[j,] ==1]
			ss <- ineighb[ineighb %in% jneighb]
			if(length(ss) == 0){
				next	
			}	
			rss <- sepsets[i,j]
			smat <- matrix(0, 2, as.integer(nchar(rss)/4))
			start <- 0
			ki <- 1
			for(k in 1:(as.integer(nchar(rss)))){
				if(substr(rss, k, k) == ","){
					smat[1,ki] <- as.numeric(substr(rss, (start+1), (k-1)))
					start <- k
				}
				if(substr(rss, k, k) == "|"){
					smat[2,ki] <- as.numeric(substr(rss, (start+1), (k-1)))
					start <- k
					ki <- ki+1
				}
	
				#ki <- (k-1)*4 + 1
				#smat[1,k] <- as.numeric(substr(rss, ki, ki))
				#smat[2,k] <- as.numeric(substr(rss, (ki+2), (ki+2)))
			}
			for(k in 1:length(ss)){
				if(!(ss[k] %in% smat)){
					ufrm2[i,ss[k]] <- 1
					ufrm2[j,ss[k]] <- 1
					ufrm2[ss[k],i] <- 0
					ufrm2[ss[k],j] <- 0		
				}	 	
			}	
			
		}
	}		
}

#r1
for(i in 1:11){
	for(j in (i+1):12){
		if(ufrm2[i,j] == ufrm2[j,i] & ufrm2[i,j] == 1){
			print(paste(i,j))
			ipar <- index[ufrm2[,i]==1]
			if(length(ipar) ==0){
				next
			}
			for(k in 1:length(ipar)){
				if(ipar[k] == j){
					next
				}
				if(ufrm2[i, ipar[k]] == 0){ #directed edge from parent to i
					ufrm2[i,j] <- 2
					ufrm2[j,i] <- 0
				}
			}
			jpar <- index[ufrm[,j] ==1]
			if(length(jpar) ==0){
				next
			}
			for(k in 1:length(jpar)){
				if(jpar[k] == i){
					next
				}
				if(ufrm2[j, jpar[k]] == 0){ #directed edge from parent to j
					ufrm2[j,i] <- 2 #had this as 1?? double check
					ufrm2[i,j] <- 0
				}
			}
		}
	}
}
#where a 2 represents -*->

#prints undirected edges
for(i in 1:11){
	for(j in (i+1):12){
		if(ufrm2[i,j] == ufrm2[j,i] & ufrm2[i,j] == 1){
			print(paste(i,j))
		}
	}
}

ufrm2s <- ufrm2
for(i in 1:12){
	for(j in 1:12){
		if(ufrm2[i,j] == 2){
			print(paste(i,j))
			ufrm2s[i,j] <- 1	
		}		
	}		
}	

g3 <- graph.adjacency(ufrm2s, mode = "directed")
plot(g3)
eg3 <- get.edgelist(g3)
E(g3)$color
for(i in 1:nrow(eg3)){
	if(ufrm2[eg3[i,1], eg3[i,2]] == 2 | ufrm2[eg3[i,2], eg3[i,1]] == 2){
		E(g3)$color[i] <- "red"	
	}
	else{
		E(g3)$color[i] <- "black"	
	}			
}	
collist <- E(g3)$color

rownames(ufrm2s) <- c("x0", "x1", "x2", "x3", "x4", "x5", "y0", "y1", "y2", "y3", "y4", "y5")
colnames(ufrm2s) <- c("x0", "x1", "x2", "x3", "x4", "x5", "y0", "y1", "y2", "y3", "y4", "y5")

g3 <- graph.adjacency(ufrm2s, mode = "directed")
E(g3)$color <- collist
plot(g3, main = paste("IC* graph of n=", nn, "and p=", prob, "and ratio of red to black = ", sum(collist=="red")/length(collist)))


g4 <- graph.adjacency(ufrm2, mode = "directed")















