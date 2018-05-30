# TODO: Add comment
# 
# Author: Paco
###############################################################################

#This is a template to do de model selection. There are libraries that you probably dont need
#I had them there and I have not cleaned those dependencies.
library(sae)
library(raster)
library(sp)
library(rgdal)

source("scripts/Functions_model_selection.R")
source("Panther_creek/Scripts/R_Scripts/Paco/Auxiliary_new.R")
source("Panther_creek/Scripts/R_Scripts/Paco/EBLUPS.R")

#read the table for the sample
Blobs_all<- read.csv("blobs_ws3_hs15/model_data.csv")

#The square of the mean elevation was useful for some variables in SWO
#I think it is necessary to make it before hand. The leaps package does not work
#make cuadratic models y~x+I(x^2)
Blobs_all$zmean_squared<-Blobs_all$zmean^2

#This predictor is a predictor Brian Wing was using and it worked well for some variables
#It coul be worth computing something similar.
#Blobs_all$Vol_Cov <-Blobs_all$Elev_mean*Blobs_all$Percentage_first_returns_above_1.50


#The script generates some diagnostic plots to help selecting a model. There are a lot of predictors
#and it is very difficult to do something completely automatic.  Other variables
#could be volumes, basal areas, or biomasses.  If you want to throw them here, it is necessary to
#add them to the Blobs table and then add their column names in the vector variables
variables <- c("log_tot_vol")

#First select the columns you will use as predictors.  I those are the one I would start with
#I think it is better not to use the intensity metrics. The intensity requires some processing
#and I doubt that is something was done in Panther_creef for both flights. Even using the
#elevation metrics and percentages or returns (plus area) I thin it is more the number of
#predictors and models is huge
predictors <- c(8:32, 36:44, 49:57)

#model_exploratory uses the leaps package functions.  Gets the (nbest) best models with 
#up to (nvmax) predictors. For the ID_SMA we could add the stand IDs. If you don't want to 
#add the stand analysis just and the blob ID and you will get some plots that are not meaningful.
#Finally, the "exps" argument is used to fit models with heteroscedaskicity, where the variance
#is proportional to the most correlated predictor in the model, 
#elevated to an exponent (in the exp vector), in the example the exponents used are:
#		0 (homocedastic model)
#		0.5(variance proportional to the most correlated predictor elevated at 0.5)
#		1(variance proportional to the most correlated predictor elevated at 0.5)
	selection_list<-c()

	for(i in variables){
		selection<-fit_candidates(Blobs_all, Blobs_all, i, predictors, nvmax = 5, 
				nbest = 5, ID_SMA = "X", force = NULL, exps = c(0,0.5,1),simplify=FALSE)
		selection_list<-c(selection_list,list(selection))
	}
	names(selection_list)<-variables
#	This is going to generate plots for all the models fitted in the previous stage
#	there will be two pdfs per dependent variable. One a plot of model summaries vs the
#	number of predictors in the model. The second pdf is the plot of the models. There are a couple things
#	to tell you in person
	lapply(selection_list,plot,output="Output/pilot_model/")
	
#	This part is filled after inspecting the  diagnostic plots of the previous stage
	candidate<- 19
	model<- 1
#This extract the selected model( models if you have more than one dependent variable in variables)
#	and generates a table with a summary of the model that gets written in the file provided in the last 
#	argument
	selected_blob_model<-tabulate_models(selection_list,candidate,model,
			"../Output/model_selection/summary_selected.txt")
