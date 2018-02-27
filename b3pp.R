
# Copyright (C) 2012-2017, Andre Falcao, University of Lisbon - LaSIGE
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
# 
# Please cite the authors in any work or product based on this material:
# 
# Martins IF, Teixeira AL, Pinheiro L, Falcao AO. 2012. A Bayesian approach to in silico blood-brain barrier penetration modeling
# J. Chem. Inf. Model., 2012, 52 (6), pp 1686-1697. DOI: 10.1021/ci300124c





load("colnames_v091x.Rout")
load("priors_v091.Rout")
load("model_v091.Rout")
sil<-capture.output(library(randomForest, quietly=T))


get_bbb<-function(params)
{
  
  b<-as.double(strsplit(params[1], " ")[[1]])	
  #create dataframe and insert data
  dat<-data.frame()
  dat<-rbind(dat,b)
  
  names(dat)<-colnames
  preds<-predict(model,dat)    #, type="prob")
  
  mw<-dat$MW
  minx<-min(priors[,1])
  maxx<-max(priors[,1])
  fac <- nrow(priors)/(maxx-minx) 
  
  #bound it for out of range Xs
  if(mw<minx) mw<-minx
  if(mw>maxx) mw<-maxx
  prior_pos <- round((mw - minx)*fac)
  prior_p<- priors[prior_pos+1, 3]   
  
  list(preds=preds, prior_p=prior_p)
}



#Here is a sample run of the model
theMol="CCO"
mol_params <-system(paste("python molattribs.py", theMol), intern=T)
res<-get_bbb(mol_params)
pred<-as.character(res$preds)
if(pred=="p"){
	cat("This Molecule should BE ABLE TO CROSS the Blood Brain Barrier")
} else{
	cat("This Molecule should NOT CROSS the Blood Brain Barrier")
}

