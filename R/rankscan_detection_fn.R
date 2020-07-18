#' @title RankScan detection
#' @description To detect and identify anomalous intervals based on rank-scans
#' @author Meng Wang
#' \email{mengw1@stanford.edu}
#' @param X A vector
#' @param alpha Level to control type I error rate (default: 0.05)
#' @param B Permutation times to get the null statistics (default: 1000)
#' @return Identified anomalous intervals under level alpha.
#' @export


#### rank scan detection
rankscan.detection.fn = function (X, alpha=0.05, B=1000) {
	n = length(X)
	Q = floor(log2(n))
	d = crt.val.rank(n, Q, alpha, B) 
	#d = 43120.513817 # n=2^15
	a = rank(X)
	names(a) = NULL
	q = floor(log2(length(X)))
	result = rank_scan_est(a, q, d)
	return(result)
}



# to get critival value for rank scan test
# only scan intervals in length 2^(1:(log2(N)-1))

GetMaxAtk <- function(X, mu.X, k){
            # k--the length of signal interval candidate	
	        N <- length(X)
	        Intv <- c(rep(1, k), rep(0, N-k))
	        return(max(convolve(X - mu.X, Intv)[1:(N-k+1)])/sqrt(k))
}

GetRankScan <- function(Rank.X, Q){
	        N <- length(Rank.X)
	        # scan intervals in dyalic lengths from 2 to [N/2]
	        #Q <- floor(log2(N))
	        max.at.k <- Q
	        for (q in 1:(Q-1)){
	        	max.at.k[q] <- GetMaxAtk(Rank.X, (N+1)/2, 2^q)
	        }
	        stats <- max(max.at.k)
            return(stats)
}


crt.val.rank = function (N, Q, alpha, B) {
	rank.scan.null <- rep(NA, B)
	for (b in 1:B){
		Rank.X.null <- sample(1:N)			
		rank.scan.null[b] <- GetRankScan(Rank.X.null, Q)
	}
	
	crt.val <- quantile(rank.scan.null, 1 - alpha)
    return(crt.val)
}



#### to get identified anomalous intervals from rank scans

rank_scan_est = function(a, q, d){
	
	# a is the rank sequence of the orignial data X
	# q is the log_2 of maximum length of the candidate intervals
	# d is the threshold for normalized sum of ranks in one interval i.e., 1/sqrt(|S|) sum_{v in S} (R_v - (N+1)/2) in the notation of the paper
	# reference for reporting the identified intervals: Jeng, X. J., Cai, T. T., and Li, H. (2010), “Optimal Sparse Segment Identification With Application in Copy Number Variation Analysis” 

	
	b = length(a)
	x = rep(0, b*(2^q+1)) # store all the intervals in length 2^(1:q)
	for (i in 1:(b-1)) {
		#print(i)
		# i is the starting index of a candidate interval
		for (j in pmin(b, i+2^(1:q)-1)) {
		    # j is its ending index
		    #print(j)
			x[(i-1)*2^q + j] = sum(a[i:j] - (b+1)/2)/sqrt(j-i+1)
			#print(x[(i-1)*2^q + j])
		}
	}
	
	
	k= which(abs(x)>d);
	i = ceiling(k/(2^q+1)); 
	j = k - (i-1)*2^q;
	list = cbind(i,j, x[k]);
	
	start = rep(0,1);
	end = rep(0,1);
	Rank_scan = rep(0,1);
	t=1;
	
	while (length(list)> 3) {
	ind = which(abs(list[,3]) == max(abs(list[,3])));
	len.ind = length(ind)
	start[t:(t+len.ind-1)] = list[ind,1];
	end[t:(t+len.ind-1)] = list[ind,2];
	Rank_scan[t:(t+len.ind-1)] = list[ind,3];
	II = c()
	for (l in 1:len.ind) {
		s = t+l-1
		II = c(II, which(list[,1]<=end[s] & list[,2]>=start[s]))
	}
	#II = which(list[,1]<=end[t] & list[,2]>=start[t]);
	list = list[-II,];
	t = t+len.ind; 
	} 
	
	if(length(list)==3) {
	start[t] = list[1];
	end[t] = list[2];
	Rank_scan[t] = list[3];
	}
	
	peaks = cbind(start, end, Rank_scan); 
	return(peaks)
}






