# TODO: Add comment
# Sep 19, 2015
# Author: maurof
###############################################################################

library(nlme)
library(leaps)
library(plyr)
library(sp)
#library(lmtest)
library(boot)
#library(rgdal)
#library(maptools)

#' Standardizes residuals
#' @export
standardized<-function(res,mcp,sigma_e,nu){
	
	res/(sigma_e*(mcp^nu))

}

#' Fits two models
#' @export
fit_two<-function(data,model_formula,weights,randomef,exp=0, model_index){
	model_formula<-formula(model_formula)
	weights<-formula(weights)
#	weights1<-do.call("nlme::varPower",list(form=weights,fixed=exp))
	if(exp==0){
		
		fixed<-try(do.call(nlme::gls,list(model=model_formula,data=data,
								method="REML",control=list(maxIter=1000,msMaxIter=1000))))
		
		mixed<-try(do.call(nlme::lme,list(fixed=model_formula,data=data,random=randomef,
								method="REML",control=list(maxIter=1000,msMaxIter=1000))))
		
	}else{
		
		fixed<-try(do.call(nlme::gls,list(model=model_formula,data=data,weights=nlme::varPower(form=weights,fixed=exp),
								method="REML",control=list(maxIter=1000,msMaxIter=1000))))
		
		mixed<-try(do.call(nlme::lme,list(fixed=model_formula,data=data,random=randomef,weights=nlme::varPower(form=weights,fixed=exp),
								method="REML",control=list(maxIter=1000,msMaxIter=1000))))
	}
	
	res<-list(fixed=fixed,mixed=mixed)
	res
	
}

#' @export
fit_four<-function(data,model_formula,weights,randomef){
	
	model_formula<-formula(model_formula)
	weights<-formula(weights)
	m00<-try(do.call(nlme::gls,list(model=formula(model_formula),data=data,
							method="REML",control=list(maxIter=1000,msMaxIter=1000))))
	
	m01<-try(do.call(nlme::lme,list(fixed=formula(model_formula),data=data,random=randomef,
							method="REML",control=list(maxIter=1000,msMaxIter=1000))))
	
	m10<-try(do.call(nlme::gls,list(model=formula(model_formula),data=data,weights=nlme::varPower(form=weights),
							method="REML",control=list(maxIter=1000,msMaxIter=1000))))
	
	m11<-try(do.call(nlme::lme,list(fixed=formula(model_formula),data=data,random=randomef,weights=nlme::varPower(form=weights),
							method="REML",control=list(maxIter=1000,msMaxIter=1000))))

	res<-list(fixed0=m00,mixed0=m01,fixed1=m10,mixed1=m11)
	res
	
}

#' @export
compare_four_models<-function(four_model_list){
	
	m00m10<-try(anova(update(four_model_list$m00,.~.,method="ML"),
					update(four_model_list$m10,.~.,method="ML"),test=TRUE))
	m00m01<-try(anova(update(four_model_list$m00,.~.,method="ML"),
					update(four_model_list$m01,.~.,method="ML"),test=TRUE))
	m10m11<-try(anova(update(four_model_list$m10,.~.,method="ML"),
					update(four_model_list$m11,.~.,method="ML"),test=TRUE))

	res<-list(m00m10=m00m10,m00m01=m00m01,m10m11=m10m11)
	
	res
}

#' @export
select_model<-function(four_model_list,mixed=TRUE){
	
	res<-compare_four_models(four_model_list)
	
	selected<-1
	
	cond_1<-(!class(res$m00m10)=="try-error")
	if(cond_1 && res$m00m10$df[1]==res$m00m10$df[2]){
		
		if(res$m00m10$AIC[2]<res$m00m10$AIC[1]){
			
			selected<-2
			
		}
		
	}else{
		
		if(cond_1 && res$m00m10$`p-value`[2]<0.05){
			
			selected<-2
			
		}
		
	}
	
	cond_2<-!cond_1 & (!class(res$m00m01)=="try-error") & mixed 
	if(cond_2 && res$m00m01$df[1]==res$m00m01$df[2]){
		
		if(res$m00m01$AIC[2]<res$m00m01$AIC[1]){
			
			selected<-3
			
		}
		
	}else{
		
		if(cond_2 && res$m00m01$`p-value`[2]<0.05){
			
			selected<-3
			
		}
		
	}
	
	cond_3<-cond_1 & (!class(res$m10m11)=="try-error") & mixed 
	if(cond_3 && res$m10m11$df[1]==res$m10m11$df[2]){
		
		if(res$m10m11$AIC[2]<res$m10m11$AIC[1]){
			
			selected<-4
			
		}
		
	}else{
		
		if(cond_3 && res$m10m11$`p-value`[2]<0.05){
			
			selected<-4
			
		}
		
	}

	four_model_list[[selected]]
}

#' @export
five_plot<-function(fitted_model,exp){
	
	if(class(fitted_model)=="try-error"){
		
		plot(-1,-1,xlim=c(0,1),ylim=c(0,1))
		text(0.5,0.5,"try-error",adj=0.5)
		plot(-1,-1,xlim=c(0,1),ylim=c(0,1))
		plot(-1,-1,xlim=c(0,1),ylim=c(0,1))
		plot(-1,-1,xlim=c(0,1),ylim=c(0,1))
		plot(-1,-1,xlim=c(0,1),ylim=c(0,1))
		return()
	}
	
	if(is.null(exp)){exp<-"NA"}
	
	if(class(fitted_model)=="lme"){
		
		plot(fitted_model$fitted[,2],residuals(fitted_model,type="pearson")
			,main=paste("exp=",exp,sep=" "),xlab="fitted values",ylab="std_residuals")
		abline(0,0,col="red")
		lines(lowess(fitted_model$fitted[,2],residuals(fitted_model,type="pearson")), col = "blue")
		
		plot(fitted_model$fitted[,2],abs(residuals(fitted_model,type="pearson")),
				main=paste("exp=",exp,sep=" "),xlab="fitted values",ylab="abs(std_residuals)")
		lines(lowess(fitted_model$fitted[,2],abs(residuals(fitted_model,type="pearson"))), col = "red")
		
		qqnorm(residuals(fitted_model,type="pearson"),main=paste("exp=",exp,sep=" "))
		qqline(residuals(fitted_model,type="pearson"),col="red")
		
		qqnorm(ranef(fitted_model)[,1],main=paste("exp=",exp,sep=" "))
		qqline(ranef(fitted_model)[,1],col="red")
		
		plot(fitted_model$fitted[,2],fitted_model$fitted[,2]+residuals(fitted_model,type="response"),
				main=paste("exp=",exp,sep=""),xlab="Predicted",ylab="Observed",
				xlim=range(c(fitted_model$fitted[,2],fitted_model$fitted[,2]+residuals(fitted_model,type="response"))),
				ylim=range(c(fitted_model$fitted[,2],fitted_model$fitted[,2]+residuals(fitted_model,type="response"))))
		abline(0,1,col="red")
		
	}else{
		
		plot(fitted_model$fitted,residuals(fitted_model,type="pearson"),
				main=paste("exp=",exp,sep=" "),xlab="fitted values",ylab="std_residuals")
		abline(0,0,col="red")
		lines(lowess(fitted_model$fitted,residuals(fitted_model,type="pearson")), col = "blue")
		
		plot(fitted_model$fitted,abs(residuals(fitted_model,type="pearson")),
				main=paste("exp=",exp,sep=" "),xlab="fitted values",ylab="abs(std_residuals)")
		lines(lowess(fitted_model$fitted,abs(residuals(fitted_model,type="pearson"))), col = "red")
		
		qqnorm(residuals(fitted_model,type="pearson"),main=paste("exp=",exp,sep=" "))
		qqline(residuals(fitted_model,type="pearson"),col="red")

		plot(-1,-1,xlim=c(0,1),ylim=c(0,1))

		plot(fitted_model$fitted,fitted_model$fitted+residuals(fitted_model,type="response"),
				main=paste("exp=",exp,sep=" "),xlab="Predicted",ylab="Observed",
				xlim=range(c(fitted_model$fitted,fitted_model$fitted+residuals(fitted_model,type="response"))),
				ylim=range(c(fitted_model$fitted,fitted_model$fitted+residuals(fitted_model,type="response"))))
		abline(0,1,col="red")
	}
		
}

#' @export
simplify_model_list<-function(x){
	
	for(i in c(1:length(x))){
		
		a<-length(x[[i]])
		to_NULL<-c()
		
		for(j in c(1:(a/2))){
			
			model_a<-try(update(x[[i]][[2*j-1]],.~.,method="ML"))
			model_b<-try(update(x[[i]][[2*j]],.~.,method="ML"))
			

			if(class(model_a)=="try-error"&class(model_b)=="try-error"){
				to_NULL<-c(to_NULL,2*j-1,2*j)
				next
				
			}else{
				
				if(class(model_a)=="try-error"){
					
					to_NULL<-c(to_NULL,2*j-1)
					next
					
				}
				if(class(model_b)=="try-error"){
					
					to_NULL<-c(to_NULL,2*j)
					next
					
				}
				
				anova_obj<-try(anova(model_a,model_b))
				
				if(anova_obj[2,'p-value']>0.05){
					
					to_NULL<-c(to_NULL,2*j)
					
				}else{
					
					to_NULL<-c(to_NULL,2*j-1)
					
				}
					
			}
			
		}
		
		to_NULL<-unique(to_NULL)
		x[[i]][to_NULL]<-NULL
		
	}
	
	x
}

#' @export
model_pars<-function(model,data_validation=NULL){
	
	if(class(model)=="try-error"){
		
		results<-data.frame(stringsAsFactors=FALSE)
		class(results)<-"try-error"
		return(results)
		
	}else{
		
		results<-data.frame(stringsAsFactors=FALSE)
		
	}
	
	variable<-all.vars(formula(model))[1]
	
	raw_residuals_training<-residuals(model,type="response")
	
	has_exp<-!is.null(model$modelStruct$varStruct)
	is_mixed<-inherits(model,"lme")
	is_spatial<-!is.null(model$modelStruct$corStruct)
	
	if(is_mixed){
		
		sigma_v<-as.numeric(nlme::VarCorr(model)["(Intercept)","StdDev"])
		coefs<-model$coefficients$fixed
		sigma_e<-model$sigma
		
	}else{
		
		sigma_v<-NA
		coefs<-model$coefficients
		sigma_e<-model$sigma
		
	}
	
	if(has_exp){
		
		if("fixed"%in%names(attributes(model$modelStruct$varStruct))){
			
			exp<-attr(model$modelStruct$varStruct,"fixed")
			
		}else{
			
			exp<-coef(model$modelStruct$varStruct,unconstrained=FALSE)
			
		}
		
		mcp<-as.character(attributes(model$modelStruct$varStruct)$formula)[2]
		
	}else{
		
		exp<-0
		mcp<-"--"
		
	}
	
	if(is_spatial){
		
		rho<-coef(model$modelStruct,unconstrained=FALSE)["corStruct.range"]
		
	}else{
		
		rho<-0
		
	}
	
	RMSE_t<-sqrt(mean(raw_residuals_training^2,na.rm=TRUE))
	Bias_t<-mean(raw_residuals_training,na.rm=TRUE)
	MAE_t<-mean(abs(raw_residuals_training),na.rm=TRUE)
	
	results[1,"variable"]<-variable
	colnames<-c()
	for(i in c(1:length(coefs))){
		
		colnames<-c(colnames,paste("beta",i-1,sep=""))
		results[1,i+1]<-coefs[i]
		
	}
	colnames(results)<-c("variable",colnames)
	results[1,"sigma_e"]<-sigma_e
	results[1,"mcp"]<-mcp
	results[1,"exp"]<-exp
	results[1,"rho"]<-rho
	results[1,"sigma_v"]<-sigma_v
	results[1,"RMSE"]<-RMSE_t
	results[1,"Bias"]<-Bias_t
	results[1,"MAE"]<-MAE_t
	results[1,"type"]<-"Training"
	results[1,"ismixed"]<-is_mixed
	results[1,"hasexp"]<-has_exp
	results[1,"isspatial"]<-is_spatial
	
	if(!is.null(data_validation)){
		
		raw_residuals_validation<-data_validation[,variable]-
				predict(model,newdata=data_validation)
		
		RMSE_v<-sqrt(mean(raw_residuals_validation^2,na.rm=TRUE))
		Bias_v<-mean(raw_residuals_validation,na.rm=TRUE)
		MAE_v<-mean(abs(raw_residuals_validation),na.rm=TRUE)
		
		results[2,"variable"]<-variable
		for(i in c(1:length(coefs))){
			
			colname<-paste("beta",i-1,sep="")
			results[2,i+1]<-coefs[i]
			
		}
		results[2,"sigma_e"]<-sigma_e
		results[2,"mcp"]<-mcp
		results[2,"exp"]<-exp
		results[2,"rho"]<-rho
		results[2,"sigma_v"]<-sigma_v
		results[2,"RMSE"]<-RMSE_v
		results[2,"Bias"]<-Bias_v
		results[2,"MAE"]<-MAE_v
		results[2,"type"]<-"Validation"
		results[2,"ismixed"]<-is_mixed
		results[2,"hasexp"]<-has_exp
		results[2,"isspatial"]<-is_spatial
		
	}
	
	results
	
}

summary.model_list<-function(x,data_validation){
	
	cont<-1
	for(i in c(1:length(x))){
		
		for(j in c(1:length(x[[i]]))){
			
			if(cont==1){
				
				res<-model_pars(x[[i]][[j]],data_validation)
				if(class(res)=="try-error"){
					
					next
					
				}
				res$candidate<-i
				res$model<-j
				cont<-cont+1
				
			}else{
				
				res2<-model_pars(x[[i]][[j]],data_validation)
				if(class(res2)=="try-error"){
					
					next
					
				}
				res2$candidate<-i
				res2$model<-j
				res<-plyr::rbind.fill(res,res2)

			}
			
		}
		
	}
	
	res
	
}

#' @export
plot.extended_leaps<-function(x,ordered=TRUE,minimize="RMSE",
		type="pdf",output="Spatial_residuals/Output/selection"){
	
	n_preds<-rownames(x$regsubsets_summary$which)
	nameY<-x$extended_summary[1,"variable"]
	hplot<-length(x$model_list[[1]])
	
	if(type=="pdf"){
		
		output1<-paste(output,nameY,".",type,sep="")
		output_comp<-paste(output,nameY,"comp.",type,sep="")
		pdf(output_comp,height=40/2.54,width=55/2.54)
		
		
	}else{
		output1<-paste(output,nameY,".",type,sep="")
		output_comp<-paste(output,nameY,"comp.",type,sep="")
		jpeg(output_comp,height=40,width=55,units="cm",res=600)
		
	}
	
	par(mfrow=c(2,2),oma=c(0,0,4,0))
	plot(as.numeric(row.names(x$regsubsets_summary$which)),x$regsubsets_summary$rsq,
			ylim=c(0,1),main="R2 vs Number of Predictors")
	plot(as.numeric(row.names(x$regsubsets_summary$which)),x$regsubsets_summary$adjr2,
			ylim=c(0,1),main="adjR2 vs Number of Predictors")
	plot(as.numeric(row.names(x$regsubsets_summary$which)),
			x$regsubsets_summary$rss,main="RSS vs Number of Predictors")
	plot(as.numeric(row.names(x$regsubsets_summary$which)),
			x$regsubsets_summary$bic,main="BIC vs Number of Predictors")
	dev.off()
	
	if(type=="pdf"){
		
		pdf(output1,height=hplot*10/2.54,width=55/2.54)
		
		
	}else{

		jpeg(output1,height=hplot*10,width=55,units="cm",res=600)
		
	}
	
	if(ordered){
		
		mins<-ddply(x$extended_summary,"candidate",summarise,
				RMSE= min(RMSE,na.rm=TRUE),
				Bias= min(Bias,na.rm=TRUE),
				MAE= min(MAE,na.rm=TRUE),
				r2= min(RMSE,na.rm=TRUE),
				adjr2= min(RMSE,na.rm=TRUE))
		
	}
	
	ordered_candidates<-order(mins[,minimize])
	par(mfrow=c(hplot,5),oma=c(0,0,4,0))
	for(i in ordered_candidates){
		cont<-1
		par(mfrow=c(hplot,5),oma=c(0,0,4,0))
		for(j in x$model_list[[i]]){
			
			
			exp<-x$extended_summary[x$extended_summary$candidate==i&
							x$extended_summary$model==cont&
							x$extended_summary$type=="Training","exp"]
			mcp<-x$extended_summary[x$extended_summary$candidate==i,"mcp"]
			mcp<-mcp[!mcp%in%c("--")][1]
			
			five_plot(j,exp)
			if(cont==1){
				
				title<-paste(as.character(formula(j)),collapse="+")
				title<-paste("candidate ",i,title," mcp=",mcp,sep="")
				print(title)
				mtext(title,line=1.5,outer=TRUE,at=0.5,adj=0.5)
				
			}
			
			cont<-cont+1
		}
		
	}
	
	dev.off()
	
}

#' @export
'[.model_list'<-function(x,i,j){
	
	x[[i]][[j]]
	
}

#' @export
'[.extended_leaps'<-function(x,i,j){

	x$model_list[i,j]
	
}

#' Fits variables
#' @export
fit_variables<-function(training,validation,variables,predictors,nvmax=8,nbest=5,ID_SMA="ID_SMA",force=NULL,
		exps=c(0:8)/4,type="jpeg",output="Spatial_residuals/Output/selection"){
	
	selection<-c()
	for(i in variables){
		
		nameY<-i

		selection_i<-fit_candidates(training,validation,nameY,predictors,nvmax=nvmax,nbest=nbest,ID_SMA=ID_SMA,force=NULL,
				exps=exps,type=type,output=output)
		
		selection<-c(selection,list(selection_i))

	}
	
	names(selection)<-variables
	
	selection<-c(selection,list(training=training,validation=validation))
	
}

#' Fits candidates
#' @export
fit_candidates<-function(training,validation,nameY,predictors,nvmax=8,nbest=5,ID_SMA="ID_SMA",force=NULL,
		exps=c(0:8)/4,simplify=TRUE){
	Y<-as.vector(training[,nameY])
	predictors<-as.matrix(training[,predictors])
	randomef<-as.formula(paste("~1|",ID_SMA,sep=""))
	candidates_obj<-leaps::regsubsets(x=predictors,y=Y,nbest=nbest,nvmax=nvmax,
			method="exhaustive",names=colnames(predictors),really.big=TRUE)
	candidates<-summary(candidates_obj)

	mcps<-c()
	cont2<-1
	models<-c()
	for(i in c(1:length(candidates$which[,1]))){
		
		print(i)
		
		model_formula<-paste(nameY,"~",paste(colnames(candidates$which[,-1])[candidates$which[i,-1]],collapse=" + "))
		vars<-c(nameY,colnames(candidates$which[,-1])[candidates$which[i,-1]])
		
		if(!is.null(force)){

			include<-!force%in%colnames(candidates$which[,-1])[candidates$which[i,-1]]
			include<-force[include]
			model_formula<-paste(c(model_formula,include),collapse=" + ")
			trainingi<-training[,c(vars,include)]
			
		}else{
			
			trainingi<-training[,vars]
			
		}
		
		cors_trainingi<-cor(trainingi)
		mcp<-rownames(cors_trainingi)[1+which(cors_trainingi[-1,1]==max(cors_trainingi[-1,1]))]
		mcps<-c(mcps,mcp)
		
		weights<-paste("~ ",mcp,sep="")
		
		if(is.null(exps)){
			
			four_models<-fit_four(training,model_formula,weights,randomef)
			models<-c(models,list(four_models))
			
		}else{
			models_i<-c()
			for(j in exps){
					
				two_models <- fit_two(training, model_formula, weights, randomef, exp=j, model_index = i)
				models_i<-c(models_i,two_models)
				
			}
			
			models<-c(models,list(models_i))
			
		}
		
	}
	
	class(models)<-c("model_list",class(models))
	if(simplify){
		
		models<-simplify_model_list(models)
		
	}
	extended_summary<-summary(models,validation)
	res<-list(training=training,validation=validation,mcps=mcps,
			regsubsets_summary=candidates,exps=exps,model_list=models,extended_summary=extended_summary)
	class(res)<-c("extended_leaps",class(res))
	res
	
}

#' @export
fit_radius<-function(selection_list,training_data,validation_data,selected_models,
		radii=7+c(11:1)*0.5,ID_SMA="ID_SMA"){
	
	variables<-names(selection_list)
	variables<-variables[-c(length(variables),length(variables)-1)]
	randomef<-as.formula(paste("~1|",ID_SMA,sep=""))
	formulas<-c()
	exps<-c()
	mcps<-c()
	
	for(i in 1:(length(selection_list)-2)){
		
		if(is.matrix(selected_models)){
			
			selected_model_i<-selected_models[i,1]
			which_model<-which(selected_model_i==selection_list[[i]]$order_candidates)
			mcp_i<-selection_list[[i]]$mcps[which_model]
			exp_i<-selected_models[i,2]
#			which_exp<-which(exp==selection_list[[i]]$models[[which_model]]$exps)
			formula_i<-formula(selection_list[[i]]$models[[which_model]]$m00)
			formulas<-c(formulas,formula_i)
			exps<-c(exps,exp_i)
			mcps<-c(mcps,mcp_i)
		
		}else{
			
			selected_model_i<-selected_models[i]
			which_model<-which(selected_model_i==selection_list[[i]]$order_candidates)
			mcp_i<-selection_list[[i]]$mcps[which_model]
			formula_i<-formula(selection_list[[i]]$models[[which_model]]$models$m00)
			formulas<-c(formulas,formula_i)
			exps<-c(exps,NULL)
			mcps<-c(mcps,mcp_i)
			
		}
		
	}
	
	list_radius<-c()
	for(radius in radii){
		
		datar<-training_data[training_data$radius==radius,]
		datavr<-validation_data[validation_data$radius==radius,]
		
		list_variables_radius<-c()
		for(i in 1:(length(selection_list)-2)){
			
			formula_i<-formulas[[i]]
			exp_i<-exps[i]
			mcp_i<-mcps[i]
			weights<-weights<-as.formula(paste("~ ",mcp_i,sep=""))
			models_var_radius<-fit_four(datar,formula_i,weights,randomef,exp_i)
			
			list_variables_radius<-c(list_variables_radius,list(models_var_radius))
			
		}
		
		names(list_variables_radius)<-variables
		list_radius<-c(list_radius,
				list(c(list_variables_radius,
						training=list(datar),validation=list(datavr))
					)
			)
		
	}
	
	names(list_radius)<-paste("radius_",radii,sep="")
	list_radius<-c(list_radius,radii=list(radii),mcps=list(mcps))
	
	list_radius
	
}

#' @export
get_selected<-function(selection_radius,plot_ID="plot_ID",ID_SMA="ID_SMA",
		mixed=FALSE,spatial=FALSE,no_aux=FALSE){
	
	radii<-selection_radius$radii
	variables<-names(selection_radius[[1]])[-c(
					length(names(selection_radius[[1]])),
					length(names(selection_radius[[1]]))-1)]
	
	mcps<-selection_radius$mcps
	
	models_radius<-list()
	for(radius in 1:length(radii)){
		
		print(radius)
		list_radius<-selection_radius[[radius]]
	
		training<-selection_radius[[radius]]$training
		validation<-selection_radius[[radius]]$validation
		
		models_var_radius<-list()
		for(variable in variables){
			
			print(variable)
			models_v_r<-list_radius[[variable]]
			
			if(no_aux){
				
				fixed_part<-as.formula(paste(variable,"~1",sep=""))
				m1<-try(do.call("gls",list(model=fixed_part,data=training,
							method="REML",control=list(maxIter=1000,msMaxIter=1000))))
				
				if(!class(m1)=="try-error"){
					
					selected_model_var_radius<-m1
					
				}
	
				if(mixed){
					
					fixed_part<-as.formula(paste(variable,"~1",sep=""))
					randomef<-as.formula(paste("~1|",ID_SMA,sep=""))
					m2<-try(do.call("lme",list(fixed=fixed_part,data=training,random=randomef,
											method="REML",control=list(maxIter=1000,msMaxIter=1000))))
					
					comp<-try(anova(m1,m2))
					
					if(!class(comp)=="try-error"){
						
						if(comp[2,6]<0.05){
						
							selected_model_var_radius<-m2
						}
					 
					}
				}
				
				

				
			}else{
				
				selected_model_var_radius<-select_model(models_v_r[1:4],mixed=mixed)
				
			}
			
			
			if(spatial){
				
				if(class(selected_model_var_radius)=="lme"){
					
					spatialcor<-as.formula(paste("~ x + y|",ID_SMA,sep=""))
					spatialcor<-corExp(form=spatialcor)
					
				}else{
					
					spatialcor<-as.formula("~ x + y ")
					spatialcor<-corExp(form=spatialcor)
					
				}
				
				model_variable_radius_sp<-try(update(selected_model_var_radius,correlation=spatialcor))
				
				comp<-try(anova(selected_model_var_radius,model_variable_radius_sp))
				
				if((!(class(comp)=="try-error"))&&comp$`p-value`[2]<0.05){
					
					selected_model_var_radius<-model_variable_radius_sp
					rm(model_variable_radius_sp)
					
				}

			}
			
			models_var_radius<-c(models_var_radius,list(selected_model_var_radius))
			
		}

		names(models_var_radius)<-variables
		models_var_radius<-c(models_var_radius,
				training=list(training),validation=list(validation))
		
		models_radius<-c(models_radius,list(models_var_radius))
		
	}
	
	names(models_radius)<-paste("radius_",radii,sep="")
	models_radius<-c(models_radius,
			radii=list(radii),variables=list(variables),
			mcps=list(mcps))
	
	models_radius
	
}

#' @export
models_pars<-function(fitted_models,
		keep_cols=c("plot_ID","MU","dir","dist","radius","x","y"),include_y=TRUE){
	
	radii<-fitted_models$radii
	variables<-fitted_models$variables
	keep_cols_a<-keep_cols
	cont<-1
	
	for(i in 1:length(radii)){
		
		radius<-radii[i]
		
		selected_models_radius<-fitted_models[[i]]
		
		for(j in 1:length(variables)){
			
			variable<-variables[j]
			
			model<-selected_models_radius[[j]]
			
			print(radius)
			print(variable)
			
			residuals_training<-residuals(model,type="pearson")
			raw_residuals_training<-residuals(model,type="response")
			datar_training<-selected_models_radius$training
			if(include_y){
				
				keep_cols<-c(keep_cols_a,variable)
				
			}
			res_training<-datar_training[,keep_cols]
			colnames(res_training)<-c(keep_cols_a,"resopnse")
			res_training<-cbind(res_training,residuals=residuals_training,
					raw_residuals=raw_residuals_training)
			res_training$variable<-variable
			
			datar_validation<-selected_models_radius$validation
			residuals_validation<-datar_validation[,variable]-
					predict(model,datar_validation)
			raw_residuals_validation<-residuals_validation
			
			resultsb<-model_pars(model,datar_training,datar_validation)
			resultsb$radius<-radius
			
			has_exp<-!is.null(model$modelStruct$varStruct)
			
			if(has_exp){
				
				exp<-coef(model$modelStruct,unconstrained=FALSE)["varStruct.power"]
				mcp<-as.character(attributes(model$modelStruct$varStruct)$formula)[2]
				sigma_e<-resultsb[1,"sigma_e"]
				residuals_validation<-standardized(residuals_validation,
						datar_validation[,mcp],sigma_e,exp)
				
			}
			
			res_validation<-datar_validation[,keep_cols]
			colnames(res_validation)<-c(keep_cols_a,"resopnse")
			res_validation<-cbind(res_validation,residuals=residuals_validation,
					raw_residuals=raw_residuals_validation)
			res_validation$variable<-variable
			
			if(cont==1){
				
				res_total_training<-res_training
				res_total_validation<-res_validation
				results<-resultsb
				
				results$radius<-radius
				
			}else{
				
				res_total_training<-rbind(res_total_training,res_training)
				res_total_validation<-rbind(res_total_validation,res_validation)
	
				results<-rbind(results,resultsb)
				
			}
			
			cont<-cont+1
		}
		
	}
	
	res<-list(results=results,
			residuals_training=res_total_training,
			residuals_validation=res_total_validation)
}

#' @export
aggregate_res<-function(data,max_radius){
	
	add1<-data[data$dist+data$radius==max_radius,]
	add1$dir<-135
	data<-rbind(data,add1)
	add1$dir<-180
	data<-rbind(data,add1)
	add1$dir<-225
	data<-rbind(data,add1)
	
	corrs_t<-data.frame(radius=numeric(),dist=numeric(),dir=numeric(),cor=numeric(),
			p.value=numeric(),significant=numeric(),LL=numeric(),UL=numeric(),variable=character())
	corrs0_t<-data.frame(radius=numeric(),dist=numeric(),dir=numeric(),cor=numeric(),
			p.value=numeric(),significant=numeric(),LL=numeric(),UL=numeric(),variable=character())
	var<-unique(data$variable)
	
	for(v in var){
		
		data_temp<-data[data$variable==v,c("plot_ID","dir","dist","radius","residuals")]
		
		corrs<-data.frame(radius=numeric(),dist=numeric(),dir=numeric(),cor=numeric(),
				p.value=numeric(),significant=numeric(),LL=numeric(),UL=numeric(),variable=character())
		corrs0<-data.frame(radius=numeric(),dist=numeric(),dir=numeric(),cor=numeric(),
				p.value=numeric(),significant=numeric(),LL=numeric(),UL=numeric(),variable=character())
		
		pred_v <- paste("pred_",v,sep="")
		res_v <- paste("res_",v,sep="")
		res_v1<-paste(res_v,"_1",sep="")
		res_v2<-paste(res_v,"_2",sep="")
		
		cont<-1
		
		for(radius in unique(data$radius)){
			
			data_tempr<-data_temp[data_temp$radius==radius,]
			max_dist<-2*max_radius-2*radius
			
			for(i in seq(from=0,to=max_dist,by=0.5)){
				
				print(v)
				print(radius)
				print(i)
				
				d1<-data_tempr[data_tempr$dist<=max_dist-i,]
				colnames(d1)<-c("plot_ID","dir","dist","radius",res_v1)
				d1$dpair<-d1$dist
				d2<-data_tempr[data_tempr$dist>=i,]
				colnames(d2)<-c("plot_ID","dir","dist","radius",res_v2)
				d2$dpair<-d2$dist-i
				d1<-merge(d1,d2,by=c("plot_ID","radius","dpair","dir"))
				rm(d2)
				
				d1<-d1[!is.na(d1[,res_v1])&!is.na(d1[,res_v2]),]
				
				res<-try(cor.test(d1[,res_v1],d1[,res_v2],method="pearson"))
				max_d_i<-max_dist-i
				res0<-try(cor.test(d1[d1$dpair==0|d1$dpair==max_d_i,res_v1],
							d1[d1$dpair==0|d1$dpair==max_d_i,res_v2],method="pearson"))
				
				if(class(res)=="try-error"){
					
					
					corrs[cont,"radius"]<-radius
					corrs[cont,"dist"]<-i
					corrs[cont,"dir"]<-361				
					corrs[cont,"cor"]<-NA
					corrs[cont,"p.value"]<-NA
					corrs[cont,"significant"]<-NA
					corrs[cont,"LL"]<-NA
					corrs[cont,"UL"]<-NA
					
					
				}else{
					
					corrs[cont,"radius"]<-radius
					corrs[cont,"dist"]<-i
					corrs[cont,"dir"]<-361				
					corrs[cont,"cor"]<-res$estimate
					corrs[cont,"p.value"]<-res$p.value
					if(is.null(res$conf.int[1])){
						
						corrs[cont,"LL"]<-NA
						
					}else{
						
						corrs[cont,"LL"]<-res$conf.int[1]
						
					}
					if(is.null(res$conf.int[2])){
						
						corrs[cont,"UL"]<-NA
						
					}else{
						
						corrs[cont,"UL"]<-res$conf.int[2]
						
					}
					
				}
				
				if(class(res0)=="try-error"){
					
					corrs0[cont,"radius"]<-radius
					corrs0[cont,"dist"]<-i
					corrs0[cont,"dir"]<-361				
					corrs0[cont,"cor"]<-NA
					corrs0[cont,"p.value"]<-NA
					corrs0[cont,"significant"]<-NA
					corrs0[cont,"LL"]<-NA
					corrs0[cont,"UL"]<-NA
					
				}else{
					
					corrs0[cont,"radius"]<-radius
					corrs0[cont,"dist"]<-i
					corrs0[cont,"dir"]<-361				
					corrs0[cont,"cor"]<-res0$estimate
					corrs0[cont,"p.value"]<-res0$p.value
					if(is.null(res0$conf.int[1])){
						
						corrs0[cont,"LL"]<-NA
						
					}else{
						
						corrs0[cont,"LL"]<-res0$conf.int[1]
						
					}
					if(is.null(res0$conf.int[2])){
						
						corrs0[cont,"UL"]<-NA
						
					}else{
						
						corrs0[cont,"UL"]<-res0$conf.int[2]
						
					}
					
				}
				
				cont<-cont+1
				
				for(j in unique(data$dir)){
					
					d3<-d1[d1$dir==j,]
					
					res<-try(cor.test(d3[,res_v1],d3[,res_v2],method="pearson"))
					res0<-try(cor.test(d3[d3$dpair==0|d3$dpair==max_d_i,res_v1],
								d3[d3$dpair==0|d3$dpair==max_d_i,res_v2],method="pearson"))
					
					if(class(res)=="try-error"){
						
						corrs[cont,"radius"]<-radius
						corrs[cont,"dist"]<-i
						corrs[cont,"dir"]<-j				
						corrs[cont,"cor"]<-NA
						corrs[cont,"p.value"]<-NA
						corrs[cont,"significant"]<-NA
						corrs[cont,"LL"]<-NA
						corrs[cont,"UL"]<-NA

						
					}else{
						
						corrs[cont,"radius"]<-radius
						corrs[cont,"dist"]<-i
						corrs[cont,"dir"]<-j				
						corrs[cont,"cor"]<-res$estimate
						corrs[cont,"p.value"]<-res$p.value
						if(is.null(res$conf.int[1])){
							
							corrs[cont,"LL"]<-NA
							
						}else{
							
							corrs[cont,"LL"]<-res$conf.int[1]
							
						}
						if(is.null(res$conf.int[2])){
							
							corrs[cont,"UL"]<-NA
							
						}else{
							
							corrs[cont,"UL"]<-res$conf.int[2]
							
						}

					}
					
					if(class(res0)=="try-error"){
						
						corrs0[cont,"radius"]<-radius
						corrs0[cont,"dist"]<-i
						corrs0[cont,"dir"]<-j				
						corrs0[cont,"cor"]<-NA
						corrs0[cont,"p.value"]<-NA
						corrs0[cont,"significant"]<-NA
						corrs0[cont,"LL"]<-NA
						corrs0[cont,"UL"]<-NA

					}else{
						
						corrs0[cont,"radius"]<-radius
						corrs0[cont,"dist"]<-i
						corrs0[cont,"dir"]<-j				
						corrs0[cont,"cor"]<-res0$estimate
						corrs0[cont,"p.value"]<-res0$p.value
						if(is.null(res0$conf.int[1])){
							
							corrs0[cont,"LL"]<-NA
							
						}else{
							
							corrs0[cont,"LL"]<-res0$conf.int[1]
							
						}
						if(is.null(res0$conf.int[2])){
							
							corrs0[cont,"UL"]<-NA
							
						}else{
							
							corrs0[cont,"UL"]<-res0$conf.int[2]
							
						}
						
					}
					
					cont<-cont+1
					
				}
				
			}
				
		}
		
		corrs[,"variable"]<-v
		corrs0[,"variable"]<-v
		
		if(v==var[1]){
			
			corrs_t<-corrs
			corrs0_t<-corrs0
			
		}else{
			
			corrs_t<-rbind(corrs_t,corrs)
			corrs0_t<-rbind(corrs0_t,corrs0)
			
		}
		
	}
	
	return(list(corrs=corrs_t,corrs0=corrs0_t))
	
}

#' @export
fit_corfunc_gls<-function(data,type="Training",start_rho=15){
	
	data<-data[!is.na(data$cor),]
	variables<-unique(data$variable)
	radii<-unique(data$radius)
	results_gls<-data.frame(variable=character(),radius=numeric(),
			rho=numeric(),type=character(),stringsAsFactors=FALSE)
	cont<-1
	for(i in variables){
		
		for(j in radii){
			
			data1<-data[data$variable==i&data$radius==j,]
			model<-try(gnls(model=cor~exp(-dist/rho), data=data1,
					start=list(rho=start_rho),
					control=list(maxIter=100,msMaxIter=100,msVerbose=TRUE)))

			results_gls[cont,"variable"]<-i
			results_gls[cont,"radius"]<-j
			results_gls[cont,"rho"]<-model$coefficients[1]
			results_gls[cont,"type"]<-type
			
			cont<-cont+1
			
		}
		
	}
	
	results_gls
}

#' @export
mixed_cor_model<-function(dist,rho,w,radius){
	
	rho1<-exp(rho)
	w1<-1/2+(atan(w)/pi)
	
	estim1<-ifelse(dist>=2*radius,0,
			(2/pi)*(acos(dist/(2*radius)))-
			(dist/(pi*radius^2))*sqrt((radius^2)-((dist/2)^2))
			)
	
	estim2<-exp((-dist/rho1))
	estim<-w1*estim1+(1-w1)*estim2
	estim
	
}

#' @export
mixed_cor_modelc<-function(dist,rho,w,radius){

	estim1<-ifelse(dist>=2*radius,0,
			(2/pi)*(acos(dist/(2*radius)))-
					(dist/(pi*radius^2))*sqrt((radius^2)-((dist/2)^2))
	)
	
	estim2<-exp((-dist/rho))
	estim<-w*estim1+(1-w)*estim2
	estim
	
}

#' @export
mixed_cor_model_gaus<-function(dist,rho,w,radius){
	
	estim1<-ifelse(dist>=2*radius,0,
			(2/pi)*(acos(dist/(2*radius)))-
					(dist/(pi*radius^2))*sqrt((radius^2)-((dist/2)^2))
	)
	
	estim2<-exp(-(dist/rho)^2)
	estim<-w*estim1+(1-w)*estim2
	estim
	
}

#' @export
grid_search<-function(start_rho,start_w,radius,max_rho,func,data){
	
	starts<-expand.grid(start_rho=start_rho,start_w=start_w)

	
	def_sol1<-1
	class(def_sol1)<-"try-error"
	def_sol2<-1
	class(def_sol2)<-"try-error"
	
	cont1<-1
	cont2<-1
	
	for(i in 1:length(starts[,1])){
		
		sol1<-try(nls(cor~func(dist,rho,w,radius=radius),
					data=data,algorithm="port",lower=c(0,0),upper=c(max_rho,1),
					weights=1/data$var_cor,
				start=list(rho=starts[i,"start_rho"],w=starts[i,"start_w"])))
		
#		print(class(sol1))
		
		sol2<-try(nls(cor~exp(-dist/rho),data=data,
				algorithm="port",lower=c(0),upper=c(max_rho),
				weights=1/data$var_cor,
				start=list(rho=starts[i,"start_rho"])))

		if(!class(sol1)=="try-error"){
			
			if(cont1==1){
				
				def_sol1<-sol1
				cont1<-cont1+1
				
			}else{
				
				
				if(deviance(def_sol1)>deviance(sol1)){
					
					def_sol1<-sol1
					
				}
				
			}	
				
		}
	
		if(!class(sol2)=="try-error"){
			
			if(cont2==1){
				
				def_sol2<-sol2
				cont2<-cont2+1
				
			}else{
				
				
				if(deviance(def_sol2)>deviance(sol2)){
					
					def_sol2<-sol2
					
				}
				
			}	
			
		}
				
	}

	comp<-try(anova(def_sol2,def_sol1))
	if(class(comp)=="try-error"){
		
		if(class(def_sol1)=="try-error"){
			
			if(class(def_sol2)=="try-error"){
				
				ret<-list(rho=NA,w=NA,range=NA)
				
			}else{
				
				sm_nls<-summary(def_sol2)
				ret<-list(rho=sm_nls$parameters["rho","Estimate"],
						w=0,
						range=3*sm_nls$parameters["rho","Estimate"])
				
			}
			
		}else{

			sm_nls<-summary(def_sol1)
			f<-function(x,rho,w,radius){
				
				0.05-func(x,rho,w,radius)
				
			}
			res<-try(uniroot(f,c(0.01,10000000000),tol=0.0001,
							radius=radius,
							rho=sm_nls$parameters["rho","Estimate"],
							w=sm_nls$parameters["w","Estimate"]))
			
			ret<-list(rho=sm_nls$parameters["rho","Estimate"],
					w=sm_nls$parameters["w","Estimate"],
					range=res$root)
				
			}
			
	}else{
		
		if(comp[2,6]<0.05){
			
			sm_nls<-summary(def_sol1)
			f<-function(x,rho,w,radius){
				
				0.05-func(x,rho,w,radius)
				
			}
			res<-try(uniroot(f,c(0.01,1000000),tol=0.0001,
							radius=radius,
							rho=sm_nls$parameters["rho","Estimate"],
							w=sm_nls$parameters["w","Estimate"]))
			
			ret<-list(rho=sm_nls$parameters["rho","Estimate"],
					w=sm_nls$parameters["w","Estimate"],
					range=res$root)
			
		}else{
			
			sm_nls<-summary(def_sol2)
			ret<-list(rho=sm_nls$parameters["rho","Estimate"],
					w=0,
					range=3*sm_nls$parameters["rho","Estimate"])
			
		}
		
	}
	
	ret<-list(ret=ret,def_sol1=def_sol1,def_sol2=def_sol2)
}

#' @export
fit_corfunc_b_gls<-function(data,type="Training",start_rho=15,start_w=0.5,criteria=1){
	
	data<-data[!is.na(data$cor),]
	variables<-unique(data$variable)
	radii<-unique(data$radius)
	results_gls<-data.frame(variable=character(),radius=numeric(),
			rho=numeric(), w=numeric(),type=character(),
			range=numeric(),sum_squares=numeric(),
			start_rho=numeric(),start_w=numeric(),stringsAsFactors=FALSE)
	cont<-1
	for(i in variables){
		
		for(j in radii){
			
			data1<-data[data$variable==i&data$radius==j,]
			radius<-j
			
			to_optim<-function(pars){
#				print(radius)
					rho1<-exp(pars[1])
					w1<-1/2+(atan(pars[2])/pi)
					
					dist<-data1$dist
					
					estim1<-ifelse(dist>=2*radius,0,
							(2/pi)*(acos(dist/(2*radius)))-
									(dist/(pi*radius^2))*sqrt((radius^2)-((dist/2)^2))
					)
					
					
					if(criteria==1){
						
						estim2<-exp((-dist/rho1))
						estim<-w1*estim1+(1-w1)*estim2
						
						result<-sum((data1$cor-estim)^2)
						
					}
					if(criteria==2){
						
						estim2<-exp((-dist/rho1))
						estim<-w1*estim1+(1-w1)*estim2
						
						result<-sum((abs(data1$cor)*(data1$cor-estim)^2))
					
						
					}
					if(criteria==3){
						
						estim2<-exp((-dist/rho1))
						estim<-w1*estim1+(1-w1)*estim2
						
						result<-sum((1/data1$var_cor)*(data1$cor-estim^2))
						
					}
					
					result
					
				}
				
			pars<-try(optim(par=c(log(start_rho),tan((start_w-0.5)*pi)),to_optim,							,
							control=list(maxit=100)))
			
			f<-function(x,rho,w,radius){
				
				0.05-mixed_cor_model(x,rho,w,radius)
				
			}
			res<-try(uniroot(f,c(0.01,1000000),tol=0.0001,
					radius=j,
					rho=pars$par[1],w=pars$par[2]))
			
			if(class(res)=="try-error"){range<-NA}else{range<-res$root}
			if(class(pars)=="try-error"){
				rho<-NA
				w<-NA
			}else{
				rho<-exp(pars$par[1])
				w<-1/2+(atan(pars$par[2])/pi)
			}
			
			results_gls[cont,"variable"]<-i
			results_gls[cont,"radius"]<-j
			results_gls[cont,"rho"]<-rho
			results_gls[cont,"w"]<-w
			results_gls[cont,"type"]<-type
			results_gls[cont,"range"]<-range
			results_gls[cont,"sum_squares"]<-pars$value
			results_gls[cont,"start_rho"]<-start_rho
			results_gls[cont,"start_w"]<-start_w
			
			cont<-cont+1
			
		}
		
	}
	
	results_gls
}

#' @export
fit_corfunc_c_gls<-function(data,type="Training",start_rho=15){
	
	data<-data[!is.na(data$cor),]
	variables<-unique(data$variable)
	radii<-unique(data$radius)
	results_gls<-data.frame(variable=character(),radius=numeric(),
			rho=numeric(), w=numeric(),type=character(),stringsAsFactors=FALSE)
	cont<-1
	for(i in variables){
		
		for(j in radii){
			
			data1<-data[data$variable==i&data$radius==j,]
			
			radius<-j
			model<-try(gnls(cor~mixed_cor_model(dist,rho,w,radius=radius),data=data1,
							start=list(rho=log(start_rho),w=tan((0.01-0.5)*pi)),
							control=list(maxIter=100,msMaxIter=100,msVerbose=TRUE)))
			results_gls[cont,"variable"]<-i
			results_gls[cont,"radius"]<-j
			results_gls[cont,"rho"]<-model$coefficients[1]
			results_gls[cont,"w"]<-model$coefficients[1]
			results_gls[cont,"type"]<-type
			
			cont<-cont+1
			
		}
		
	}
	
	results_gls
}

#' @export
fit_corfunc_d_gls<-function(start_rho=15,start_w,max_rho,func=mixed_cor_modelc,data=data,type="Residuals"){
	
	data<-data[!is.na(data$cor),]
	variables<-unique(data$variable)
	radii<-unique(data$radius)
	results_gls<-data.frame(variable=character(),radius=numeric(),
			rho=numeric(), w=numeric(),type=character(),range=numeric(),
			stringsAsFactors=FALSE)
	list_results<-c()
	cont<-1
	for(i in variables){
		
		for(j in radii){
			
			data1<-data[data$variable==i&data$radius==j,]
			
			radius<-j
			sol<-grid_search(start_rho,start_w,radius=radius,max_rho=max_rho,func,data=data1)
			sol<-list(variable=i,radius=j,sol=sol)
			list_results<-c(list_results,sol=list(sol))
			results_gls[cont,"variable"]<-i
			results_gls[cont,"radius"]<-j
			results_gls[cont,"rho"]<-sol$sol$ret$rho
			results_gls[cont,"w"]<-sol$sol$ret$w
			results_gls[cont,"range"]<-sol$sol$ret$range
			results_gls[cont,"type"]<-type
			
			cont<-cont+1
			
		}
		
	}
	
	results_gls<-list(results_gls=results_gls,list_results=list_results)
}

#Other functions
#' @export
add_levels<-function(data,ID_SMA,all_levels){
	
	data[,ID_SMA]<-as.factor(data[,ID_SMA])
	data[,paste(ID_SMA,"_old",sep="")]<-data[,ID_SMA]
	newlevels<-all_levels[!all_levels%in%levels(data[,ID_SMA])]	
	levels(data[,ID_SMA])<-c(levels(data[,ID_SMA]),newlevels)
	data
	
	
}

#' @export
plot_errores<-function(lme.obj,variable,density,up=TRUE,left=TRUE,type="p",a=1.2,b=1.2,c=1.2,d=1.5){
	
#	if(up){par(oma=c(2,0,4,0))}else{par(oma=c(4,0,2,0))}
	par(mar=c(5,5,5,0.5),cex.axis=1.4)
	qqnorm(lme.obj$coefficients$random[[1]],
			main="",pch=20,ylab="",xlab="")
	qqline(lme.obj$coefficients$random[[1]],col="red")
	mtext("Q-Q plot stand level\nrandom effects",3,line=0.5,cex=a,outer=FALSE)
	mtext("Theoretical Quantiles",1,line=3,cex=a,outer=FALSE)
	mtext("Sample Quantiles",2,line=3,cex=a,outer=FALSE)
	box("figure")
	par(mar=c(5,0.5,5,5))
	qqnorm(residuals(lme.obj,type=type),
			main="",
			yaxt="n",pch=20,cex.main=b,ylab="",xlab="")
	axis(side=4)
	mtext("Q-Q plot unit level\nrandom effects",3,line=0.5,cex=a,outer=FALSE)
	mtext("Theoretical Quantiles",1,line=3,cex=a,outer=FALSE)
	mtext("Sample Quantiles",4,line=3,cex=a,outer=FALSE)
	qqline(residuals(lme.obj,type=type),col="red")
	box("figure")
	
	if(type=="p"){
		
		res_tit<-"Standarized Residuals"
		
	}else{
		
		
		res_tit<-"Residuals"
		
	}
	
	if(left){
		
		par(mar=c(5,5,5,1))
		plot(lme.obj$fitted[,2],residuals(lme.obj,type=type),main="",
				cex.main=b,xlab="",ylab="",pch=20,yaxt="n",ylim=c(-3,3))
		axis(side=2)
		mtext("Predicted Value",1,line=3,cex=a,outer=FALSE)
		abline(h=0,col="red")
		mtext(res_tit,2,line=3,cex=a,outer=FALSE)
		box("figure")
		
	}else{
		
		par(mar=c(5,1,5,5))
		plot(lme.obj$fitted[,2],residuals(lme.obj,type=type),main="",
				cex.main=b,xlab="",ylab="",pch=20,yaxt="n",ylim=c(-3.5,3.5))
		axis(side=4)
		mtext("Predicted Value",1,line=3,cex=a,outer=FALSE)
		abline(h=0,col="red")
		
		mtext(res_tit,4,line=3,cex=a,outer=FALSE)
		box("figure")
		
	}
	
	mtext("Predicted vs Residuals",3,line=3,cex=a,outer=FALSE)
	
	if(up){
		
		if(left){at<-0.25}else{at<-0.75}
		
		
		mtext(paste("Model for ",variable,sep=""),3,outer=TRUE,at=at,adj=0.5,cex=d)
		
		
	}else{
		
		if(left){at<-0.25}else{at<-0.75}
		
		mtext(paste("Model for ",variable,sep=""),3,outer=TRUE,at=at,adj=0.5,cex=d)
		
	}
	
}

#' @export
plot_comb<-function(m1,m2,variable1,variable2,title,output,width,height,res=res,units="cm",type="p"){
	
#	postscript(output,height=height,width=width)
	layout(matrix(c(1,2,4,5,3,3,6,6),ncol=4,byrow=TRUE))
	par(oma=c(4,5,4,4))
	for(i in c(1:2)){
		
#		b<-unlist(a[i],recursive=FALSE)
#		
#		lme.obj<-b$lme
#		Dens<-b$densidad
#		Variable<-b$Variable
#		preds<-b$Preds
		
		
		if(i==1){
			
			plot_errores(m1,variable=variable1,density="",up=TRUE,left=TRUE,type=type)
			
		}
		
		if(i==2){
			
			plot_errores(m2,variable=variable2,density="",up=TRUE,left=FALSE,type=type)
			
		}
		
#		if(i==3){
#			
#			plot_errores(lme.obj,Dens,up=FALSE,left=TRUE)
#			
#		}
#		
#		if(i==4){
#			
#			plot_errores(lme.obj,Dens,up=FALSE,left=FALSE)
#			
#		}
		
	}
#	mtext(title,3,outer=TRUE,line=3,adj=0.5)
#	box("outer")
#	dev.off()
	
}

#selection of models
#' @export
prepare_model<-function(selection,variables,selected_model_number,ID_SMA="MU",
		write.df=NULL,tol=1e-6,opt="optim"){
	
	{
		selection<-selection[variables]
		selected_model_number<-selected_model_number[variables]
		
		variables<-names(selection)
		selected_model_number <- selected_model_number[order(variables)]
		variables<-variables[order(variables)]
		
		sel_exps<-c()
		varcovs<-c()
		names(selected_model_number)<-variables
		random_formula<-paste("~variable-1|",ID_SMA,spe="")
		random_formula<-as.formula(random_formula)
		corr_formula<-paste("~ varindex|",ID_SMA,"/rowindex")
		corr_formula<-as.formula(corr_formula)
		
		Plots<-selection[[1]]$data
		n_plots<-length(Plots[,1])
		
		varindex<-c()
		cont<-1
		res_var_exps<-c()
		res_var_covariate<-c()
		res_model_number<-c()
		var_exps<-c()
		var_covariates<-c()
		var_covariates_names<-c()
		rows<-c()
		variabs<-c()
		RMSEs_univariate<-c()
		univariate<-c()
		for(i in variables){
			
			modelnumber<-which(selection[[i]]$models_stats$model_cont==selected_model_number[i])
			res_model_number<-c(res_model_number,modelnumber)
			linear_model<-selection[[i]]$mixed_models[[modelnumber]]
			univariate<-c(univariate,list(linear_model))
			RMSEs_univariate<-c(RMSEs_univariate,sqrt(mean(linear_model$residuals[,2]^2,na.rm=TRUE)))
			var_exp<-c(selection[[i]]$models_stats[modelnumber,"exp"])
			res_var_exps<-c(res_var_exps,var_exp)
			var_covariate<-selection[[i]]$models_stats[modelnumber,"varformula"]
			var_covariate<-strsplit(var_covariate,"~ ")[[1]][2]
			res_var_covariate<-c(res_var_covariate,var_covariate)
			var_exps<-c(var_exps,rep(var_exp,length.out=length(Plots[,1])))
			var_covariates_names<-c(var_covariates_names,
					rep(var_covariate,length.out=length(Plots[,1])))
			var_covariates<-c(var_covariates,Plots[,var_covariate]^var_exp)
			rows<-c(rows,1:length(Plots[,1]))
			variabs<-c(variabs,rep(i,length.out=length(Plots[,1])))
			
			varindex<-c(varindex,cont)
			if(cont==1){
				
				vars<-colnames(linear_model$data)
				rand_effects<-data.frame(unique(levels(linear_model$data[,ID_SMA])[linear_model$data[,ID_SMA]]),
						linear_model$coefficients$random[[1]])
				resdf<-data.frame(linear_model$residuals[,2])
				cont<-cont+1
				
			}else{
				
				newvars<-colnames(linear_model$data)
				vars<-c(vars,newvars[which(!newvars%in%vars)])
				rand_effects<-cbind(rand_effects,
						data.frame(linear_model$coefficients$random[[1]]))
				resdf<-cbind(resdf,
						data.frame(linear_model$residuals[,2]))
				cont<-cont+1
			}
			
			
			
		}
		
		colnames(rand_effects)<-c(ID_SMA,variables)
		initmatcor_rand<-cor(rand_effects[,-1])
		colnames(resdf)<-variables
		initcor_res<-cor(resdf)
		
		names(univariate)<-variables
		initcor<-c()
		for(i in 1:length(variables)-1){
			
			if((i+1)<length(variables)){
				
				for(j in (i+1):length(variables)){
					
					
					initcor<-c(initcor,initcor_res[i,j])
					
				}
				
				
			}else{
				
				initcor<-c(initcor,initcor_res[i,i+1])
				
			}
			
			
		}
		
		tomerge<-data.frame(variable=variabs,rowindex=rows,var_exps=var_exps,
				var_covariate=var_covariates,var_covariates_name=var_covariates_names)
		
		
		names(RMSEs_univariate)<-variables
		
		vars1<-which(vars%in%variables)
		vars2<-which(!vars%in%variables)
		vars1name<-vars[which(vars%in%variables)]
		vars2name<-vars[which(!vars%in%variables)]
		vars3name<-vars2name[-which(vars2name==ID_SMA)]
		vars<-c(vars,ID_SMA)
		
		multivariatedf<-as.matrix(cbind(rep(1,times=length(Plots[,ID_SMA])),Plots[,vars3name]))
		expand<-diag(length(variables))
		colnamesX<-paste(paste("V_",1:length(variables),"_",sep=""),
				rep(c("Intercept",vars3name),each=length(variables)),sep="")
		multivariatedf<-as.data.frame(multivariatedf%x%expand)
		colnames(multivariatedf)<-colnamesX
		
		cont<-1
		fixed_formula<-"valuey~-1+variable"
		pick_beta<-c()
		for(j in vars3name){
			
			for(i in variables){
				
				modelnumber<-which(selection[[i]]$models_stats$model_cont==selected_model_number[i])
				linear_model<-selection[[i]]$mixed_models[[modelnumber]]
				preds<-colnames(linear_model$data)
				
				if(j%in%preds){
					
					fixed_formula<-paste(fixed_formula,"+",colnamesX[cont+length(variables)])
					pick_beta<-c(pick_beta,colnamesX[cont+length(variables)])
				}
				cont<-cont+1
				
			}
			
			
			
			
		}
		
		
		
		fixed_formula<-as.formula(fixed_formula)
		multivariatedf[,ID_SMA]<-as.factor(rep(Plots[,ID_SMA],each=length(variables)))
		
		cont<-1
		for(i in variables){
			
			which0<-rep(0,length(variables))
			which0[cont]<-1
			
			
			if(cont==1){
				
				y<-Plots[,i]
				y<-rep(y,each=length(variables))
				y<-y*which0
				cont<-cont+1
				
			}else{
				
				y2<-Plots[,i]
				y2<-rep(y2,each=length(variables))
				y2<-y2*which0
				y<-y+y2
				cont<-cont+1
				
				
			}
			
			
			
		}
		multivariatedf[,"valuey"]<-y
		multivariatedf[,"variable"]<-rep(variables,times=length(Plots[,ID_SMA]))
		multivariatedf[,"rowindex"]<-rep(c(1:length(Plots[,ID_SMA])),each=length(variables))
		
		multivariatedf<-merge(multivariatedf,tomerge,by=c("variable","rowindex"))
		multivariatedf$rowindex<-as.factor(multivariatedf$rowindex)
		
		
		multivariatedf<-multivariatedf[order(levels(multivariatedf[,ID_SMA])[multivariatedf[,ID_SMA]],
						multivariatedf$rowindex,
						multivariatedf$variable),]
		
		multivariatedf[,"varindex"]<-rep(c(1:length(variables)),times=length(Plots[,ID_SMA]))
		
		varindexSMA<-paste("varindex",ID_SMA,sep="")
		rowindexSMA<-paste("rowindex",ID_SMA,sep="")
		rowindexvarindexSMA<-paste("varindex_","rowindex",ID_SMA,sep="")
		
		multivariatedf[,varindexSMA]<-interaction(multivariatedf[,"varindex"],
				multivariatedf[,ID_SMA])
		multivariatedf[,varindexSMA]<-as.factor(multivariatedf[,varindexSMA])
		
		multivariatedf[,rowindexSMA]<-interaction(multivariatedf[,"rowindex"],
				multivariatedf[,ID_SMA])
		multivariatedf[,rowindexSMA]<-as.factor(multivariatedf[,rowindexSMA])
		
		
		multivariatedf$rowindexvariable<-with(multivariatedf,interaction(rowindex,variable))
		multivariatedf$rowindexvariable<-as.factor(multivariatedf$rowindexvariable)
		
		multivariatedf[,rowindexvarindexSMA]<-interaction(multivariatedf[,"rowindex"],
				multivariatedf[,"varindex"],multivariatedf[,ID_SMA])
		multivariatedf[,rowindexvarindexSMA]<-as.factor(multivariatedf[,rowindexvarindexSMA])
		
		multivariatedf$var_covariates_name<-as.factor(multivariatedf$var_covariates_name)
		multivariatedf$variable<-as.factor(multivariatedf$variable)
		
#	multivariatedf<-multivariatedf[order(levels(multivariatedf[,ID_SMA])[multivariatedf[,ID_SMA]],
#					multivariatedf$rowindex,
#					multivariatedf$variable),]
		
	}
	
	
	multivariate_cor<-multivariate<-lme(fixed_formula,
			data=multivariatedf,
			random=random_formula,
			correlation=corSymm(value =initcor,form =corr_formula),
			weights=varComb(varIdent(value=1,~1|variable),nlme::varPower(form=~var_covariate,fixed=1)),
			control=list(opt=opt,maxIter=1000,
					msMaxIter=1000,
					niterEM=2000,
					msMaxEval=1000,
					msMaxIter=1000,
					msVerbose=TRUE,
					tolerance=tol))
	
	multivariate_no_cor<-multivariate<-lme(fixed_formula,
			data=multivariatedf,
			random=random_formula,
			weights=varComb(varIdent(value=1,~1|variable),nlme::varPower(form=~var_covariate,fixed=1)),
			control=list(opt=opt,maxIter=1000,
					msMaxIter=1000,
					niterEM=2000,
					msMaxEval=1000,
					msMaxIter=1000,
					msVerbose=TRUE,
					tolerance=tol))
	
	getcols<-c(paste("V_",1:length(variables),"_Intercept",sep=""),pick_beta)
	X<-as.matrix(multivariatedf[,getcols])
	pick_beta<-c(paste("variable",variables,sep=""),pick_beta)
	beta<-multivariate_cor$coefficients$fixed[pick_beta]
	
	
	nrows_Z<-dim(multivariatedf)[1]
	ncols_Z<-length(variables)*length(levels(multivariatedf[,ID_SMA]))
	dim_R<-ncols_Z
	
	Z<-matrix(0,ncol=ncols_Z,nrow=nrows_Z)
	
	
	sigma<-multivariate_cor$sigma
	sigmas<-as.vector(sigma*attr(multivariate_cor$modelStruct$varStruct$A,"weights"))*
			as.vector(attr(multivariate_cor$modelStruct$varStruct$B,"covariate"))
	
	R1<-as.matrix(corMatrix(multivariate_cor$modelStruct$corStruct)[[1]])
	R<-diag(1,n_plots)%x%R1
	R<-R*(sigmas%*%t(sigmas))
	R2<-NA
#	R2<-sigmas%*%t(sigmas)
	
	for(i in 1:nrows_Z){
		
		area<-multivariatedf[i,ID_SMA]
		area<-which(area==levels(multivariatedf[,ID_SMA]))
		variable1<-multivariatedf[i,"variable"]
		variable1<-which(variable1==variables)
		nvars<-length(variables)
		Z[i,(area-1)*nvars+variable1]<-1
		
#		for(j in 1:nrows_Z){
#			
#			variable2<-multivariatedf[i,"variable"]
#			variable2<-which(variable2==variables)
#			R2[i,j]<-sigmas[i]*sigmas[j]*R1[variable1,variable2]
#			
#			
#		}
		
	}
#	R2<-(diag(1,n_plots)%x%R1)*R2
	G<-getVarCov(multivariate_cor)
	G<-G[,paste("variable",variables,sep="")]
	G<-G[paste("variable",variables,sep=""),]
	G<-diag(length(levels(multivariatedf[, ID_SMA])))%x%G
	
	V<-R+Z%*%G%*%t(Z)
	
	V_inv<-solve(V)
	
	W<-solve(solve(G)+t(Z)%*%solve(R)%*%Z)
	
	y<-multivariatedf[,"valuey"]
	
	residuals<-data.frame(residuals=multivariate_cor$residuals[,2])
	multivariatedf<-cbind(multivariatedf,residuals)
	
	RMSEs_multivariate<-c()
	for(i in variables){
		
		RMSEs_multivariate<-c(RMSEs_multivariate,
				sqrt(mean(multivariatedf[multivariatedf$variable==i,"residuals"]^2,na.rm=TRUE)))
		
		
	}
	
	names(RMSEs_multivariate)<-variables
	RMSEs<-data.frame(RMSEs_univariate=RMSEs_univariate,RMSEs_multivariate=RMSEs_multivariate)
	rownames(RMSEs)<-variables
	
	
	res<-list(variables=variables,multivariate_cor=multivariate_cor,
			multivariate_no_cor=multivariate_no_cor,
			univariate=univariate,
			multivariatedf=multivariatedf,n_plots=n_plots,
			y=y,X=X,Z=Z,beta=beta,G=G,R=R,R1=R1,R2=R2,V=V,V_inv=V_inv,W=W,
			fixed_formula=fixed_formula,
			random_formula=random_formula,
			corr_formula=corr_formula,
			initcor=initcor,ID_SMA=ID_SMA,
			beta=beta,getcols=getcols,variables=variables,
			selected_model_number=selected_model_number,
			rand_effects=rand_effects,cor_rand=initmatcor_rand,
			residuals_univariate=resdf,cor_res=initcor_res,
			RMSEs=RMSEs,
			predictors=vars3name,
			res_var_exps=res_var_exps,
			res_var_covariate=res_var_covariate,
			res_model_number=res_model_number
	)
	
	
	
	
}

#' @export
prepare_model_optim<-function(selection,selected_model_number,ID_SMA="MU",
		write.df=NULL,tol=1e-6){
	
	{
		
#	variables<-names(selection)
		selected_model_number <- selected_model_number[order(variables)]
		variables<-variables[order(variables)]
		names(selected_model_number)<-variables
		
		sel_exps<-c()
		varcovs<-c()
		names(selected_model_number)<-variables
		random_formula<-paste("~variable-1|",ID_SMA,spe="")
		random_formula<-as.formula(random_formula)
		corr_formula<-paste("~ varindex|",ID_SMA,"/rowindex")
		corr_formula<-as.formula(corr_formula)
		
		Plots<-selection[[1]]$data
		n_plots<-length(Plots[,1])
		
		varindex<-c()
		cont<-1
		var_exps<-c()
		var_covariates<-c()
		var_covariates_names<-c()
		rows<-c()
		variabs<-c()
		RMSEs_univariate<-c()
		univariate<-c()
		init_sigma_v<-c()
		init_sigma_e<-c()
		for(i in variables){
			
			modelnumber<-which(selection[[i]]$models_stats$model_cont==selected_model_number[i])
			
			linear_model<-selection[[i]]$mixed_models[[modelnumber]]
			init_sigma_e<-c(init_sigma_e,linear_model$sigma)
			init_sigma_v<-c(init_sigma_v,getVarCov(linear_model)[1,1])
			univariate<-c(univariate,list(linear_model))
			RMSEs_univariate<-c(RMSEs_univariate,sqrt(mean(linear_model$residuals[,2]^2,na.rm=TRUE)))
			var_exp<-c(selection[[i]]$models_stats[modelnumber,"exp"])
			var_covariate<-selection[[i]]$models_stats[modelnumber,"varformula"]
			var_covariate<-strsplit(var_covariate,"~ ")[[1]][2]
			var_exps<-c(var_exps,rep(var_exp,length.out=length(Plots[,1])))
			var_covariates_names<-c(var_covariates_names,
					rep(var_covariate,length.out=length(Plots[,1])))
			var_covariates<-c(var_covariates,Plots[,var_covariate]^var_exp)
			rows<-c(rows,1:length(Plots[,1]))
			variabs<-c(variabs,rep(i,length.out=length(Plots[,1])))
			
			varindex<-c(varindex,cont)
			if(cont==1){
				
				vars<-colnames(linear_model$data)
				rand_effects<-data.frame(unique(levels(linear_model$data[,ID_SMA])[linear_model$data[,ID_SMA]]),
						linear_model$coefficients$random[[1]])
				resdf<-data.frame(linear_model$residuals[,2])
				cont<-cont+1
				
			}else{
				
				newvars<-colnames(linear_model$data)
				vars<-c(vars,newvars[which(!newvars%in%vars)])
				rand_effects<-cbind(rand_effects,
						data.frame(linear_model$coefficients$random[[1]]))
				resdf<-cbind(resdf,
						data.frame(linear_model$residuals[,2]))
				cont<-cont+1
			}
			
			
			
		}
		
		colnames(rand_effects)<-c(ID_SMA,variables)
		initmatcor_v<-cor(rand_effects[,-1])
		colnames(resdf)<-variables
		initmatcor_e<-cor(resdf)
		
		names(univariate)<-variables
		initcor_v<-c()
		initcor_e<-c()
		for(i in 1:length(variables)-1){
			
			if((i+1)<length(variables)){
				
				for(j in (i+1):length(variables)){
					
					
					initcor_e<-c(initcor_e,initmatcor_e[i,j])
					initcor_v<-c(initcor_v,initmatcor_v[i,j])
				}
				
				
			}else{
				
				initcor_e<-c(initcor_e,initmatcor_e[i,i+1])
				initcor_v<-c(initcor_v,initmatcor_v[i,j])
				
			}
			
			
		}
		
		tomerge<-data.frame(variable=variabs,rowindex=rows,var_exps=var_exps,
				var_covariate=var_covariates,var_covariates_name=var_covariates_names)
		
		
		names(RMSEs_univariate)<-variables
		
		vars1<-which(vars%in%variables)
		vars2<-which(!vars%in%variables)
		vars1name<-vars[which(vars%in%variables)]
		vars2name<-vars[which(!vars%in%variables)]
		vars3name<-vars2name[-which(vars2name==ID_SMA)]
		vars<-c(vars,ID_SMA)
		
		multivariatedf<-as.matrix(cbind(rep(1,times=length(Plots[,ID_SMA])),Plots[,vars3name]))
		expand<-diag(length(variables))
		colnamesX<-paste(paste("V_",1:length(variables),"_",sep=""),
				rep(c("Intercept",vars3name),each=length(variables)),sep="")
		multivariatedf<-as.data.frame(multivariatedf%x%expand)
		colnames(multivariatedf)<-colnamesX
		
		cont<-1
		fixed_formula<-"value~-1+variable"
		pick_beta<-c()
		for(j in vars3name){
			
			for(i in variables){
				
				modelnumber<-which(selection[[i]]$models_stats$model_cont==selected_model_number[i])
				linear_model<-selection[[i]]$mixed_models[[modelnumber]]
				preds<-colnames(linear_model$data)
				
				if(j%in%preds){
					
					fixed_formula<-paste(fixed_formula,"+",colnamesX[cont+length(variables)])
					pick_beta<-c(pick_beta,colnamesX[cont+length(variables)])
				}
				cont<-cont+1
				
			}
			
			
			
			
		}
		
		
		
		fixed_formula<-as.formula(fixed_formula)
		multivariatedf[,ID_SMA]<-as.factor(rep(Plots[,ID_SMA],each=length(variables)))
		
		cont<-1
		for(i in variables){
			
			which0<-rep(0,length(variables))
			which0[cont]<-1
			
			
			if(cont==1){
				
				y<-Plots[,i]
				y<-rep(y,each=length(variables))
				y<-y*which0
				cont<-cont+1
				
			}else{
				
				y2<-Plots[,i]
				y2<-rep(y2,each=length(variables))
				y2<-y2*which0
				y<-y+y2
				cont<-cont+1
				
				
			}
			
			
			
		}
		multivariatedf[,"value"]<-y
		multivariatedf[,"variable"]<-rep(variables,times=length(Plots[,ID_SMA]))
		multivariatedf[,"rowindex"]<-rep(c(1:length(Plots[,ID_SMA])),each=length(variables))
		
		multivariatedf<-merge(multivariatedf,tomerge,by=c("variable","rowindex"))
		multivariatedf$rowindex<-as.factor(multivariatedf$rowindex)
		
		multivariatedf<-multivariatedf[order(levels(multivariatedf[,ID_SMA])[multivariatedf[,ID_SMA]],
						multivariatedf$rowindex,
						multivariatedf$variable),]
		
		multivariatedf[,"varindex"]<-rep(c(1:length(variables)),times=length(Plots[,ID_SMA]))
		
		varindexSMA<-paste("varindex",ID_SMA,sep="")
		rowindexSMA<-paste("rowindex",ID_SMA,sep="")
		rowindexvarindexSMA<-paste("varindex_","rowindex",ID_SMA,sep="")
		
		multivariatedf[,varindexSMA]<-interaction(multivariatedf[,"varindex"],
				multivariatedf[,ID_SMA])
		multivariatedf[,varindexSMA]<-as.factor(multivariatedf[,varindexSMA])
		
		multivariatedf[,rowindexSMA]<-interaction(multivariatedf[,"rowindex"],
				multivariatedf[,ID_SMA])
		multivariatedf[,rowindexSMA]<-as.factor(multivariatedf[,rowindexSMA])
		
		
		multivariatedf$rowindexvariable<-with(multivariatedf,interaction(rowindex,variable))
		multivariatedf$rowindexvariable<-as.factor(multivariatedf$rowindexvariable)
		
		multivariatedf[,rowindexvarindexSMA]<-interaction(multivariatedf[,"rowindex"],
				multivariatedf[,"varindex"],multivariatedf[,ID_SMA])
		multivariatedf[,rowindexvarindexSMA]<-as.factor(multivariatedf[,rowindexvarindexSMA])
		
		multivariatedf$var_covariates_name<-as.factor(multivariatedf$var_covariates_name)
		multivariatedf$variable<-as.factor(multivariatedf$variable)
		
		
		
		
		getcols<-c(paste(paste("V_",1:length(variables),"_",sep=""),
						"Intercept",sep=""),pick_beta)
		X<-as.matrix(multivariatedf[,getcols])
		pick_beta<-c(paste("variable",variables,sep=""),pick_beta)
#	beta<-multivariate$coefficients$fixed[pick_beta]
		
		nrows_Z<-dim(multivariatedf)[1]
		ncols_Z<-length(variables)*length(levels(multivariatedf[,ID_SMA]))
		dim_R<-ncols_Z
		
		Z<-matrix(0,ncol=ncols_Z,nrow=nrows_Z)
		for(i in 1:nrows_Z){
			
			area<-multivariatedf[i,ID_SMA]
			area<-which(area==levels(multivariatedf[,ID_SMA]))
			variable1<-multivariatedf[i,"variable"]
			variable1<-which(variable1==variables)
			nvars<-length(variables)
			Z[i,(area-1)*nvars+variable1]<-1
			
#		for(j in 1:nrows_Z){
#			
#			variable2<-multivariatedf[i,"variable"]
#			variable2<-which(variable2==variables)
#			R2[i,j]<-sigmas[i]*sigmas[j]*R1[variable1,variable2]
#			
#			
#		}
			
		}
		
		var_covariates<-multivariatedf[,"var_covariate"]
		
		
		initpars<-c(init_sigma_v,init_sigma_e,initcor_v,initcor_e)
		lower<-rep(0,2*length(variables))
		upper<-rep(1e+18,2*length(variables))
		lower<-c(lower,rep(-1,(length(variables)*(length(variables)-1))))
		upper<-c(upper,rep(1,(length(variables)*(length(variables)-1))))
		
	}
	
	to_optim<-function(pars,X,Z,y,var_covariates,variables,n_plots){
		
		sigmas_v<-sqrt(pars[1:length(variables)])
		sigmas_e<-sqrt(pars[(length(variables)+1):(2*length(variables))])
		cors<-pars[(2*length(variables)+1):length(pars)]
		cors_v<-cors[1:(length(cors)/2)]
		cors_e<-cors[-c(1:(length(cors)/2))]
		
		
		R1<-matrix(0,nrow=length(variables),ncol=length(variables))
		G1<-R1
		cont<-1
		for(i in 1:length(variables)){
			
			for(j in i:length(variables)){
				
				if(i==j){
					
					R1[i,j]<-sigmas_e[i]*sigmas_e[j]
					G1[i,j]<-sigmas_v[i]*sigmas_v[j]
				}else{
					
					
					R1[i,j]<-sigmas_e[i]*sigmas_e[j]*cors_e[cont]
					R1[j,i]<-sigmas_e[i]*sigmas_e[j]*cors_e[cont]
					
					G1[i,j]<-sigmas_v[i]*sigmas_v[j]*cors_v[cont]
					G1[j,i]<-sigmas_v[i]*sigmas_v[j]*cors_v[cont]
					
					cont<-cont+1
					
				}
				
				
			}
			
		}
		
		
		R<-as.matrix(var_covariates%*%t(var_covariates))
		R_cors<-diag(1,n_plots)%x%R1
		R<-R*R_cors
		
		G<-diag(length(levels(multivariatedf[, ID_SMA])))%x%G1
		
		V<-try(as.spam(R+Z%*%G%*%t(Z)))
		
		V_inv<-try(solve(V))
		XtVX_inv<-try(solve(t(X)%*%V_inv%*%X))
#		print(determinant(as.spam(XtVX_inv))$modulus)
		P<-try(as.spam(V_inv-(V_inv%*%X%*%XtVX_inv%*%t(X)%*%V_inv)))
		det_V_Spam<-try(determinant(V)$modulus)
		
		if(any(c(class(V),class(V_inv),class(XtVX_inv),class(P),class(det_V_Spam))=="try-error")){
			
			ll<-1e+18
			return(ll)
			
		}
		
		
		
		if(any(c(sigmas_e,sigmas_v)<=0)){
			
			ll<-1e+18
			return(ll)
		}else{
			
			ll<-det_V_Spam[1]-t(y)%*%P%*%y
			ll<--ll[1,1]
			
		}
		
		
		if(any(abs(c(cors_v,cors_e))>1)){
			
			ll<-1e+18
			return(ll)
		}
		
		ll
	}
	
	sol_LBFGSB<-optim(par=initpars,to_optim,X=X,Z=Z,y=y,
			var_covariates=var_covariates,
			variables=variables,n_plots=n_plots,
			method="L-BFGS-B",upper=upper,lower=lower,
			control=list(trace=6))
	
	
	sol_SANN<-optim(par=initpars,to_optim,X=X,Z=Z,y=y,
			var_covariates=var_covariates,
			variables=variables,n_plots=n_plots,
			method="SANN",
			control=list(trace=6))
	
	sol_NM<-optim(par=initpars,to_optim,X=X,Z=Z,y=y,
			var_covariates=var_covariates,
			variables=variables,n_plots=n_plots,
			method="Nelder-Mead",
			control=list(trace=6))
	
	
	
	
	
	
	
	
	
#	R2<-sigmas%*%t(sigmas)
	
	
#	R2<-(diag(1,n_plots)%x%R1)*R2
	
	
	W<-solve(solve(G)+t(Z)%*%solve(R)%*%Z)
	
	residuals<-data.frame(residuals=multivariate$residuals[,2])
	multivariatedf<-cbind(multivariatedf,residuals)
	
	RMSEs_multivariate<-c()
	for(i in variables){
		
		RMSEs_multivariate<-c(RMSEs_multivariate,
				sqrt(mean(multivariatedf[multivariatedf$variable==i,"residuals"]^2,na.rm=TRUE)))
		
		
	}
	
	names(RMSEs_multivariate)<-variables
	RMSEs<-data.frame(RMSEs_univariate=RMSEs_univariate,RMSEs_multivariate=RMSEs_multivariate)
	rownames(RMSEs)<-variables
	
	res<-list(variables=variables,multivariate=multivariate,univariate=univariate,
			multivariatedf=multivariatedf,
			X=X,Z=Z,beta=beta,G=G,R=R,R1=R1,R2=R2,V=V,V_inv=V_inv,
			fixed_formula=fixed_formula,
			random_formula=random_formula,
			corr_formula=corr_formula,
			initcor=initcor,ID_SMA=ID_SMA,
			beta=beta,variables=variables,
			selected_model_number=selected_model_number,
			rand_effects=rand_effects,residuals_univariate=resdf,
			RMSEs=RMSEs,
			predictors=vars3name	
	)
	
	
	
	
}

#' @export
create_predictions_multivariate<-function(res_multivariate=res_multivariate,data=data){	
	
	ID_SMA<-res_multivariate$ID_SMA
	predictors<-res_multivariate$predictors
	variables<-res_multivariate$variables
	var_covariates<-res_multivariate$res_var_covariate
	multivariatedf<-res_multivariate$multivariatedf
	multivariate_cor<-res_multivariate$multivariate_cor
	
	y<-res_multivariate$y
	X<-res_multivariate$X
	Z<-res_multivariate$Z
	G<-res_multivariate$G
	V_inv<-res_multivariate$V_inv
	W<-res_multivariate$W
	
	
	data[,"Id_row"]<-c(1:length(data[,1]))+res_multivariate$n_plots
	data[,"Intercept"]<-rep(1,length.out=length(data[,1]))
	data<-data[order(data[,"Id_row"]),]
	
	
	IDs_SMA<-rep(data[,ID_SMA],each=length(variables))
	L<-as.matrix(data[,c("Intercept",predictors)])
	new_names<-paste(paste("V_",1:length(variables),"_",sep=""),
			rep(colnames(L),each=length(variables)),sep="")
	new_names[1:length(variables)]<-paste("variable",variables,sep="")
	L<-L%x%diag(1,length(variables))
	colnames(L)<-new_names
	rownames(L)<-rep(data[,"Id_row"],each=length(variables))
	beta<-res_multivariate$beta
	L<-L[,names(beta)]
	
	var_covs<-c()
	for(i in 1:length(variables)){
		
		var_covs<-c(var_covs,data[,var_covariates[i]])
		
	}
	
	var_covs<-data.frame(variable=rep(variables,each=length(data[,1])),
			var_covs=var_covs,
			Id_row=rep(data[,"Id_row"],times=length(variables)))
	
	var_covs<-var_covs[order(var_covs[,"Id_row"],var_covs[,"variable"]),]
	
	var_covs<-var_covs[,"var_covs"]
	
	L_df<-data.frame(Id_row=rep(data[,"Id_row"],each=length(variables)),
			variable=rep(variables,times=length(data[,"Id_row"])),L)
	
	L_df<-merge(L_df,data[,c("Id_row",ID_SMA)],by="Id_row")
	
	nrows_M<-dim(L)[1]
	ncols_M<-length(variables)*length(levels(multivariatedf[,ID_SMA]))
	M<-matrix(0,ncol=ncols_M,nrow=nrows_M)
	for(i in 1:nrows_M){
		
		area<-L_df[i,ID_SMA]
		area<-which(area==levels(multivariatedf[,ID_SMA]))
		variable1<-L_df[i,"variable"]
		variable1<-which(variable1==variables)
		nvars<-length(variables)
		M[i,(area-1)*nvars+variable1]<-1
		
		
	}
	
	sigma<-multivariate_cor$sigma
	
	sigmas<-as.vector(sigma*
							attr(multivariate_cor$modelStruct$varStruct$A,"weights")[1:length(variables)][variables])*
			var_covs
	
	sigmas<-data.frame(Id_row=rep(data[,"Id_row"],times=length(variables)),sigmas=sigmas)
	
	Q1<-as.matrix(corMatrix(multivariate_cor$modelStruct$corStruct)[[1]])
	G4s<-by(sigmas,INDICES=list(Id_row=as.factor(sigmas[,"Id_row"])),function(x){
				
				
				ret<-x$sigmas%*%t(x$sigmas)*Q1
				
			},simplify=FALSE)
	
	list_ID_SMA<-list(IDs_SMA)
	names(list_ID_SMA)<-ID_SMA
	Q_ID_SMA<-aggregate(L[,1],by=list_ID_SMA[1],function(x){
				
				res<-sum (x)
				
			})
	
	K<-M%*%G%*%t(Z)%*%V_inv
	var_cov_beta<-solve(t(X)%*%V_inv%*%X)
	
	preds<-L%*%beta
	
	preds<-preds+K%*%(y-X%*%beta)
	
	preds<-data.frame(preds=preds[,1],
			Id_row=rep(data[,"Id_row"],each=length(variables)),
			variable=rep(variables,times=length(data[,1])))
	
	maps<-by(preds,INDICES=list(variable=preds[,"variable"]),function(x){
				
				return(x[,c("preds","Id_row")])
				
			})
	
	
	
	data[,"Area"]<-1
	univariates<-c()
	cont<-1
	for(i in variables){
		
		map<-maps[[i]]
		colnames(map)<-c(paste("preds_",i,sep=""),"Id_row")
		
#		univariate<-eblups(res_multivariate$univariate[[i]],data,"Id_row")
#		univariates<-c(univariates,list(univariate))
		data<-merge(data,map,by="Id_row",all.x=TRUE)
		
		
	}
	
	names(univariates)<-variables
	
	
	
	M_df<-data.frame(Id_row=rep(data[,"Id_row"],each=length(variables)),
			variable=rep(variables,times=length(data[,"Id_row"])),M)
	
	M_df<-merge(M_df,data[,c("Id_row",ID_SMA)],by="Id_row")
	
	List_ID_SMA_vars<-list(name=L_df[,ID_SMA],name2=L_df[,"variable"])
	names(List_ID_SMA_vars)<-c(ID_SMA,"variable")
	
	L_g<-aggregate(L_df[,-c(1,2,dim(L_df)[2])],by=List_ID_SMA_vars,mean,na.rm = TRUE)
	L_g<-L_g[order(L_g[,ID_SMA],L_g[,"variable"]),]
	
	M_g<-aggregate(M_df[,-c(1,2,dim(M_df)[2])],by=List_ID_SMA_vars,mean,na.rm = TRUE)
	M_g<-M_g[order(M_g[,ID_SMA],M_g[,"variable"]),]
	
	cols<-c(ID_SMA,paste("preds_",variables,sep=""))
	
	List_ID_SMAs<-list(name=data[,ID_SMA])
	names(LIst_ID_SMAs)<-ID_SMA
	predictionsSMAs<-aggregate(data[,cols],by=List_ID_SMAs,mean,na.rm = TRUE)
	
	
	for(i in variables){
		
		tomerge_pixels<-univariates[[i]]$eblups$pixels[,c(3:9)]
		tomerge_pixels[,"Ids_2"]<-as.numeric(
				levels(tomerge_pixels[,"Ids"])[tomerge_pixels[,"Ids"]])
		names(tomerge_pixels[,c(1:5,7)])<-paste("univariate",i,colnames(tomerge_pixels[,c(1:5,7)]),sep="_")
		data<-merge(data,tomerge_pixels,by.x=c("Id_row"),by.y=c("Ids_2"),all.x=TRUE)
		
		tomerge_g<-univariates[[i]]$eblups$stands[,c(1,3:7,9)]
		names(tomerge_g[,c(2:7)])<-paste("univariate",i,colnames(tomerge_pixels[,c(2:7)]),sep="_")
		predictionsSMAs<-merge(predictionsSMAs,tomerge_g,by.x=c(ID_SMA),by.y=c(ID_SMA),all.x=TRUE)
		
		
	}
	
	indexes_G1s<-expand.grid(V1=c(1:length(variables)),V2=c(1:length(variables)))
	indexes_G2s<-paste("G2",indexes_G1s[,1],indexes_G1s[,2],sep="_")
	indexes_G4s<-paste("G4",indexes_G1s[,1],indexes_G1s[,2],sep="_")
	indexes_G1s<-paste("G1",indexes_G1s[,1],indexes_G1s[,2],sep="_")
	
	data[,indexes_G1s]<-0
	data[,indexes_G2s]<-0
	data[,indexes_G4s]<-0
	
	predictionsSMAs[,indexes_G1s]<-0
	predictionsSMAs[,indexes_G2s]<-0
	
	
	
	for(i in predictionsSMAs[,ID_SMA]){
		
		M_i<-M_g[M_g[,ID_SMA]==i,]
		M_i<-M_i[order(M_i[,ID_SMA],M_i[,"variable"]),]
		M_i<-as.matrix(M_i[,-c(1,2)])
		
		L_i<-L_g[L_g[,ID_SMA]==i,]
		L_i<-L_i[order(L_i[,ID_SMA],L_i[,"variable"]),]
		L_i<-as.matrix(L_i[,-c(1,2)])
		
		
		
		G1_g<-M_i%*%W%*%t(M_i)
		G2_g<-(L_i-M_i%*%G%*%t(Z)%*%V_inv%*%X)
		G2_g<-G2_g%*%var_cov_beta%*%t(G2_g)
		
		for(j in c(1:length(variables))){
			
			for(k in c(1:length(variables))){
				
				col_G1<-paste("G1",j,k,sep="_")
				col_G2<-paste("G2",j,k,sep="_")
				
				predictionsSMAs[predictionsSMAs[,ID_SMA]==i,col_G1]<-G1_g[j,k]
				predictionsSMAs[predictionsSMAs[,ID_SMA]==i,col_G2]<-G2_g[j,k]
				
			}
			
			
		}
		
		
		
	}
	
	save.image("F:/Paco/Corvallis/FIA_meeting/Presentation/multivariate_inside_presentation.RData")
	write.table(predictionsSMAs,"F:/Paco/Corvallis/FIA_meeting/Presentation/SMALLareas_presentation.txt",row.names=FALSE)
	
	for(i in data[,"Id_row"]){
		
		M_i<-M_df[M_df[,"Id_row"]==i,]
		M_i<-M_i[order(M_i[,"Id_row"],M_i[,"variable"]),]
		M_i<-as.matrix(M_i[,-c(1,2,dim(M_i)[2])])
		
		L_i<-L_df[L_df[,"Id_row"]==i,]
		L_i<-L_i[order(L_i[,"Id_row"],L_i[,"variable"]),]
		L_i<-as.matrix(L_i[,-c(1,2,dim(L_i)[2])])
		
		print(i)
		
		G1_g<-M_i%*%W%*%t(M_i)
		G2_g<-(L_i-M_i%*%G%*%t(Z)%*%V_inv%*%X)
		G2_g<-G2_g%*%var_cov_beta%*%t(G2_g)
		
#		G4_g<-G4s[[i]]
		
		for(j in c(1:length(variables))){
			
			for(k in c(1:length(variables))){
				
				col_G1<-paste("G1",j,k,sep="_")
				col_G2<-paste("G2",j,k,sep="_")
				col_G4<-paste("G4",j,k,sep="_")
				
				data[data[,"Id_row"]==i,col_G1]<-G1_g[j,k]
				data[data[,"Id_row"]==i,col_G2]<-G2_g[j,k]
				data[data[,"Id_row"]==i,col_G4]<-G4_g[j,k]
				
			}
			
			
		}
		
		
		
	}
	
	save.image("F:/Paco/Corvallis/FIA_meeting/Presentation/multivariate_inside.RData")
	
	res<-list(L=L,M=M,Q_ID_SMA=Q_ID_SMA,pixels=data,stands=predictionsSMAs,univariates=univariates)
	res
	
}

#' @export
paco_pdf <- function(model_frame, nvmax = 5, nbest = 5, ID_SMA) {
	for(i in variables){
  		selection <- fit_candidates(model_frame, model_frame, i, predictors, nvmax = 5, 
                            nbest = 5, ID_SMA = "OI_KEY", force = NULL, exps = c(0,0.5,1),simplify=FALSE)
  		selection_list <- c(selection_list,list(selection))
	}

	names(selection_list) <- variables
	return(selection_list)
}