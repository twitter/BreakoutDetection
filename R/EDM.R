breakout = function(Z, min.size = 30, method = 'amoc', ...){

	#function used to scale observations to the interval [0,1]
	f = function(x) (x-min(x))/(max(x)-min(x))
	
	if(class(Z)%in%c('numeric','integer') || ncol(Z) == 1)
		Zcounts = f(Z)
	else if(!('timestamp'%in%names(Z) && 'count'%in%names(Z)))
		stop("The supplied data must either be a vector, or a data.frame which has columns named 'timestamp' and 'count'.")
	else
		Zcounts = f(Z$count)
	
	#capture additional passed arguments
	argList = list(...)


	#Argument checking
	if(min.size < 2)
		stop("min.size must be an integer greater than 2.")
	if(!(method %in% c('amoc','multi')))
		stop("method must be either 'amoc' or 'multi'.")
	
	#For single change detection check for necessary arguments
	#if they aren't present set default value
	if(method == 'amoc'){
		multi = F
		val = argList[['alpha']]
		if(is.null(val) || val > 2 || val <= 0)
			alpha = 2
		else
			alpha = val
		
		val = argList[['sig.lvl']]
		if(is.null(val) || val <= 0 || val >= 1)
			sig.lvl = 0.05
		else
			sig.lvl = val
		
		val = argList[['nperm']]
		if(is.null(val) || val < 0)
			nperm = 0
		else
			nperm = floor(val)
		
		val = argList[['exact']]
		if(is.null(val) || !is.logical(val))
			exact = T
		else
			exact = val
	}
	#For multiple change detection check for the degree argument
	#For beta and percent check is performed later
	else if(method == 'multi'){
		multi = T
		degree = argList[['degree']]
		if(!(degree %in% c(0,1,2)) || is.null(degree))
			degree = 1

		beta = argList[['beta']]
		percent = argList[['percent']]
	}
	#Select the correct EDM subroutine to run
	if(!multi && exact)
		Analysis = EDMX

	else if(!multi && !exact)
		Analysis = EDM_tail

	#Check for penalty type (beta or percent)
	else if(multi && ( !is.null(beta) || is.null(percent) )){
		if(is.null(beta))
			beta = 0.008
		penalty = beta
		Analysis = EDM_multi
	}
	
	else if(multi && is.null(beta) && !is.null(percent)){
		penalty = percent
		Analysis = EDM_percent
	}

	#Perform analysis. Report runnint time and change locations
	if(!multi){
		p1 = proc.time()
		retList = Analysis(Zcounts, min.size, alpha)
		if(nperm == 0){
			p2 = proc.time()
			retList$time = as.numeric((p2-p1)[3])
			retList$pval = NA
		}
		else{
			over = 1
			for(i in 1:nperm){
				Zcounts = sample(Zcounts)
				stat = Analysis(Zcounts, min.size, alpha)
				if(stat$stat >= retList$stat)
					over = over + 1
			}
			p2 = proc.time()
			retList$time = as.numeric((p2-p1)[3])
			retList$pval = over/(nperm+1)
			if(retList$pval > sig.lvl)
				retList$loc = NULL
		}
	}

	else{
		p1 = proc.time()
		retList = Analysis(Zcounts, min.size, penalty, degree)
		p2 = proc.time()
		retList$time = as.numeric((p2-p1)[3])
		retList$pval = NA
	}

	#Check to see if plotting was desired and generate plot if necessary
	plotting = argList[['plot']]
	title = argList[['title']]
	xlab = argList[['xlab']]
	ylab = argList[['ylab']]

	if(!is.null(plotting) && plotting == T){
		dateTime = T
		if(class(Z)%in%c('numeric','integer') || ncol(Z) == 1){
			dateTime = F
			Z = data.frame(timestamp=1:length(Z), count = Z)
		}

		g = ggplot2::ggplot(Z, ggplot2::aes(x=timestamp, y=count)) + ggplot2::theme_bw() + 
		    ggplot2::theme(panel.grid.minor=ggplot2::element_blank(), panel.grid.major=ggplot2::element_blank())
		
		if(is.null(xlab))
			xlab = " "
		if(is.null(ylab))
			ylab = " "
		if(is.null(title))
			title = " "
			
		g = g + ggplot2::xlab(xlab) + ggplot2::ylab(ylab) + ggplot2::ggtitle(title)

		g = g + ggplot2::geom_line()
		if(!is.null(retList$loc)&& length(retList$loc)>0){
			v = retList$loc
			v = c(0,v)
			for(j in 2:length(v)){
				M = mean(Z$count[(v[j-1]+1):v[j]])
				df2 = data.frame(Z$timestamp[v[j]], Z$timestamp[v[j]], -Inf, M)
				names(df2) = c('x','xend','y','yend')
				g = g + ggplot2::geom_segment(data=df2,ggplot2::aes(x=x,y=y,xend=xend,yend=yend,color='2'),linetype=2,size=1.2)
				g = g + ggplot2::guides(color=FALSE)
			}
		}
		if(dateTime)
			g = g + scale_x_datetime(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
		else
			g = g + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))

		retList$plot = g
	}
	
	return(retList)
}
