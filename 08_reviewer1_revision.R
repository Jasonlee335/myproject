# Reviewer-requested analyses; Steps 1-7 and their outputs are preserved.
options(stringsAsFactors = FALSE, scipen = 999)
stopifnot(all(vapply(c("survey", "survival"), requireNamespace, logical(1), quietly = TRUE)))
suppressPackageStartupMessages(library(survey))
suppressPackageStartupMessages(library(survival))
options(survey.lonely.psu = "adjust")
set.seed(20261028)
project <- Sys.getenv("DXA_PROJECT_DIR", file.path(getwd(), "DXA_transportability_results"))
out <- file.path(project, "reviewer1_revision", "analysis_results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
B <- as.integer(Sys.getenv("REVISION_BOOTSTRAP_REPS", "1000"))
stopifnot(is.finite(B), B >= 100)
wmean <- function(x,w) {
  ok <- is.finite(x) & is.finite(w) & w>0
  if (!any(ok)) return(NA_real_)
  sum(x[ok]*w[ok])/sum(w[ok])
}
ci <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)<.9*B) return(c(NA_real_,NA_real_))
  unname(quantile(x,c(.025,.975),type=6))
}
capacity_rule <- function(s,w,capacity=.2) {
  ok <- is.finite(s) & is.finite(w) & w>0
  o <- which(ok)[order(s[ok],decreasing=TRUE)]
  threshold <- s[o[which(cumsum(w[o])>=capacity*sum(w[o]))[1]]]
  above <- sum(w[ok & s>threshold]); at <- sum(w[ok & s==threshold])
  fraction <- min(1,max(0,(capacity*sum(w[ok])-above)/at))
  list(threshold=threshold,fraction=fraction,
       referral=ifelse(s>threshold,1,ifelse(s==threshold,fraction,0)))
}
apply_rule <- function(s,rule) ifelse(s>rule$threshold,1,ifelse(s==rule$threshold,rule$fraction,0))
metrics <- function(y,r,w) {
  ok <- y %in% c(0,1) & is.finite(r) & is.finite(w) & w>0
  y<-y[ok];r<-r[ok];w<-w[ok]
  tp<-sum(w*y*r);cases<-sum(w*y);selected<-sum(w*r)
  c(referral_share=selected/sum(w),sensitivity=tp/cases,ppv=tp/selected,
    examinations_per_case=selected/tp,prevalence=cases/sum(w))
}
weighted_auc <- function(y,s,w) {
  ok<-y %in% c(0,1) & is.finite(s) & is.finite(w) & w>0
  if(length(unique(y[ok]))!=2) return(NA_real_)
  a<-aggregate(cbind(case_w=w[ok]*y[ok],control_w=w[ok]*(1-y[ok])),list(score=s[ok]),sum)
  a<-a[order(a$score),]
  sum(a$case_w*(cumsum(a$control_w)-.5*a$control_w))/(sum(a$case_w)*sum(a$control_w))
}
rep_weights <- function(d) {
  des<-svydesign(ids=~psu,strata=~strata,weights=~survey_weight,nest=TRUE,data=d)
  as.matrix(weights(as.svrepdesign(des,type="bootstrap",replicates=B,mse=TRUE),type="analysis"))
}
write_out <- function(x,name) write.csv(x,file.path(out,name),row.names=FALSE,na="")
im<-readRDS(file.path(project,"data_clean","imaging_cross_country_analysis_ready.rds"))
periods<-c(development="NHANES 2005-2010",temporal_validation_1="NHANES 2013-2014",
  temporal_validation_2="NHANES 2017-Mar 2020",historical_geographic_validation="KNHANES 2008-2011",
  contemporary_temporal_validation="KNHANES 2024")
im$period<-unname(periods[im$analysis_role])
im$age_band<-cut(im$age,c(50,60,70,80),right=FALSE,labels=c("50-59","60-69","70-79"))
im$age5<-cut(im$age,seq(50,80,5),right=FALSE)
im$risk<--im$osta_raw
sets<-split(im,factor(im$period,levels=unname(periods)))
rw<-lapply(sets,rep_weights)
source_name<-unname(periods["development"]);src<-sets[[source_name]]
source_rule<-capacity_rule(src$risk,src$survey_weight)
source_reps<-lapply(seq_len(B),function(b) capacity_rule(src$risk,rw[[source_name]][,b]))
q<-tapply(src$survey_weight,src$age5,sum);q<-q/sum(q)
age_rows<-aggregate_rows<-drift_rows<-std_rows<-list()
idx<-ai<-di<-si<-0L
for(pn in names(sets)) {
  message("Age/capacity/bootstrap: ",pn)
  d<-sets[[pn]];w<-d$survey_weight;wr<-rw[[pn]]
  r_age<-numeric(nrow(d))
  for(g in levels(d$age_band)) {
    k<-d$age_band==g;r_age[k]<-capacity_rule(d$risk[k],w[k])$referral
  }
  strategies<-list(published=as.numeric(d$osta<=-1),source20=apply_rule(d$risk,source_rule),
    pooled20=capacity_rule(d$risk,w)$referral,age_stratified20=r_age)
  for(ref in c("common_US","local")) {
    y<-if(ref=="common_US") d$osteoporosis_fn_us else d$osteoporosis_fn_local
    for(sn in names(strategies)) {
      r<-strategies[[sn]]
      bs<-matrix(NA_real_,B,5,dimnames=list(NULL,names(metrics(y,r,w))))
      ba<-array(NA_real_,c(B,3,4))
      for(b in seq_len(B)) {
        wb<-wr[,b]
        rb<-switch(sn,published=r,
          source20=apply_rule(d$risk,if(pn==source_name) source_reps[[b]] else source_rule),
          pooled20=capacity_rule(d$risk,wb)$referral,
          age_stratified20={
            rr<-numeric(nrow(d))
            for(g in levels(d$age_band)) {
              k<-d$age_band==g;rr[k]<-capacity_rule(d$risk[k],wb[k])$referral
            }
            rr
          })
        bs[b,]<-metrics(y,rb,wb)
        for(gg in seq_along(levels(d$age_band))) {
          k<-d$age_band==levels(d$age_band)[gg]
          mm<-metrics(y[k],rb[k],wb[k])
          ba[b,gg,]<-c(mm["referral_share"],mm["sensitivity"],
            sum(wb[k]*rb[k])/sum(wb*rb),weighted_auc(y[k],d$risk[k],wb[k]))
        }
      }
      point<-metrics(y,r,w)
      for(m in names(point)) {
        ai<-ai+1L;interval<-ci(bs[,m])
        aggregate_rows[[ai]]<-data.frame(period=pn,reference=ref,strategy=sn,metric=m,
          estimate=unname(point[m]),lower=interval[1],upper=interval[2],n=nrow(d),cases=sum(y==1))
      }
      for(gg in seq_along(levels(d$age_band))) {
        g<-levels(d$age_band)[gg];k<-d$age_band==g;mm<-metrics(y[k],r[k],w[k])
        values<-c(referral_share=mm[["referral_share"]],sensitivity=mm[["sensitivity"]],
          selected_composition=sum(w[k]*r[k])/sum(w*r),auc=weighted_auc(y[k],d$risk[k],w[k]))
        for(m in seq_along(values)) {
          idx<-idx+1L;interval<-ci(ba[,gg,m])
          age_rows[[idx]]<-data.frame(period=pn,reference=ref,strategy=sn,age_band=g,
            metric=names(values)[m],estimate=values[m],lower=interval[1],upper=interval[2],
            n=sum(k),cases=sum(y[k]==1),population_share=sum(w[k])/sum(w),
            case_share=sum(w[k]*y[k])/sum(w*y),sparse=sum(y[k]==1)<20)
        }
      }
    }
  }
  # Direct standardization conditions on a fixed source age distribution and cutoff.
  for(sn in c("published","source20")) {
    r<-strategies[[sn]]
    by5<-vapply(levels(d$age5),function(g) wmean(r[d$age5==g],w[d$age5==g]),numeric(1))
    std<-sum(q*by5)
    bs<-vapply(seq_len(B),function(b) {
      by5<-vapply(levels(d$age5),function(g) wmean(r[d$age5==g],wr[d$age5==g,b]),numeric(1))
      sum(q*by5)
    },numeric(1))
    interval<-ci(bs);si<-si+1L
    std_rows[[si]]<-data.frame(period=pn,strategy=sn,observed_referral=wmean(r,w),
      source_age_standardized_referral=std,lower=interval[1],upper=interval[2],
      note="Six 5-year groups; source age distribution and cutoff fixed")
  }
  if(pn!=source_name) {
    cond<-joint<-matrix(NA_real_,B,2,dimnames=list(NULL,c("referral_share","sensitivity")))
    for(b in seq_len(B)) {
      cond[b,]<-metrics(d$osteoporosis_fn_us,strategies$source20,wr[,b])[1:2]
      joint[b,]<-metrics(d$osteoporosis_fn_us,apply_rule(d$risk,source_reps[[b]]),wr[,b])[1:2]
    }
    for(mode in c("target_only_conditional","source_and_target_joint")) {
      mat<-if(mode=="target_only_conditional") cond else joint
      for(m in colnames(mat)) {
        interval<-ci(mat[,m]);di<-di+1L
        drift_rows[[di]]<-data.frame(period=pn,bootstrap=mode,metric=m,
          estimate=metrics(d$osteoporosis_fn_us,strategies$source20,w)[m],
          lower=interval[1],upper=interval[2],successful_replicates=sum(is.finite(mat[,m])))
      }
    }
  }
}
ages<-do.call(rbind,age_rows);agg<-do.call(rbind,aggregate_rows)
write_out(ages,"R1_age_stratified_performance.csv")
write_out(agg,"R1_capacity_tradeoffs.csv")
write_out(do.call(rbind,std_rows),"R1_age_standardized_referral.csv")
write_out(do.call(rbind,drift_rows),"R1_joint_source_target_bootstrap.csv")
write_out(data.frame(source_threshold_osta=-source_rule$threshold,tie_fraction=source_rule$fraction,
                     bootstrap_replicates=B),"R1_source_threshold.csv")

ch<-readRDS(file.path(project,"data_clean","charls_analysis_ready.rds"))
ch$risk_z<-as.numeric(scale(-ch$osta_raw))
o<-order(ch$osta_raw,as.character(ch$id));ch$high20<-0L
ch$high20[o[seq_len(round(.2*nrow(ch)))]]<-1L
ch<-ch[complete.cases(ch[,c("followup_years","incident_hip","risk_z","high20","baseline_fall",
  "rural","stroke","arthritis","current_smoking","ever_drinking","community_id","bmi")]),]
factors<-"baseline_fall + rural + stroke + arthritis"
forms<-c(clinical=paste("Surv(followup_years,incident_hip) ~",factors),
  clinical_plus_osta=paste("Surv(followup_years,incident_hip) ~ risk_z +",factors),
  clinical_plus_age_weight=paste("Surv(followup_years,incident_hip) ~ age + weight_kg +",factors),
  clinical_plus_age_weight_splines=paste("Surv(followup_years,incident_hip) ~ splines::ns(age,3) + splines::ns(weight_kg,3) +",factors),
  clinical_plus_high20=paste("Surv(followup_years,incident_hip) ~ high20 +",factors))
fits<-lapply(forms,function(f) coxph(as.formula(f),data=ch,cluster=community_id,x=TRUE))
conc<-do.call(rbind,lapply(names(fits),function(nm) {
  f<-fits[[nm]]
  data.frame(model=nm,n=f$n,events=f$nevent,apparent_harrell_c=unname(summary(f)$concordance[1]),
    predictors=forms[nm],note="In-sample descriptive C-index; not external validation")
}))
write_out(conc,"R1_charls_age_weight_comparators.csv")
primary<-do.call(rbind,lapply(c("clinical_plus_osta","clinical_plus_high20"),function(nm) {
  f<-fits[[nm]];term<-if(nm=="clinical_plus_osta") "risk_z" else "high20"
  s<-summary(f)$coefficients[term,];beta<-unname(s["coef"]);se<-unname(s["robust se"])
  data.frame(exposure=term,hr=exp(beta),lower=exp(beta-1.96*se),upper=exp(beta+1.96*se),
    p_value=2*pnorm(-abs(beta/se)),n=f$n,events=f$nevent)
}))
write_out(primary,"R1_charls_primary_exact_p.csv")
saveRDS(list(age_results=ages,capacity=agg,models=fits,bootstrap_reps=B),file.path(out,"R1_analysis_objects.rds"))
writeLines(capture.output(sessionInfo()),file.path(out,"R1_sessionInfo.txt"))
inputs<-file.path(project,"data_clean",c("imaging_cross_country_analysis_ready.rds","charls_analysis_ready.rds"))
write_out(data.frame(file=basename(inputs),md5=unname(tools::md5sum(inputs))),"R1_input_manifest.csv")
cat("Reviewer 1 analyses completed with",B,"bootstrap replicates.\nResults:",out,"\n")
