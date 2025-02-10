#多种模型的批量单因素cox回归分析
colnames(adjustedtm)<-gsub("-",".",colnames(adjustedtm))
comsam<-intersect(colnames(adjustedtm),rownames(allcox))
adjustedtm<-adjustedtm[,comsam]
#-----------------------------b-----------------------
rownames(b)<-b$gene
comsam<-intersect(rownames(adjustedtm),rownames(b))
b<-cbind(adjustedtm[comsam,],b)
b<-b[,-3047]
write.xlsx(b,'b.xlsx')
scores <- colSums(b[, 1:ncol(b)-1] * b$coef)
b<-rbind(b,scores)
rownames(b)[10]<-'b'
b<-b[10,]
b<-b[,-3047]
b<-t(b)
save(b,file='b.RData')

#------------------------------------d----------------------
rownames(d)<-d$gene
comsam<-intersect(rownames(adjustedtm),rownames(d))
d<-cbind(adjustedtm[comsam,],d[,1])
write.xlsx(d,'d.xlsx')
scores <- colSums(d[, 1:ncol(d)-1] * d$coef)
d<-rbind(d,scores)
rownames(d)[6]<-'d'
d<-d[6,]
d<-d[,-3047]
d<-t(d)
save(d,file='d.RData')

#-------------------------e--------------------
rownames(e)<-e$gene
comsam<-intersect(rownames(adjustedtm),rownames(e))
e<-cbind(adjustedtm[comsam,],e[,1])
write.xlsx(e,'e.xlsx')
scores <- colSums(e[, 1:ncol(e)-1] * e$coef)
e<-e[,-3047]
e<-rbind(e,scores)
rownames(e)[6]<-'e'
e<-e[6,]
e<-t(e)
save(e,file='e.RData')

#------------------------f----------------------
rownames(f)<-f$gene
comsam<-intersect(rownames(adjustedtm),rownames(f))
f<-cbind(adjustedtm[comsam,],f[,1])
write.xlsx(f,'f.xlsx')
scores <- colSums(f[, 1:ncol(f)-1] * f$coef)
f<-f[,-3047]
f<-rbind(f,scores)
rownames(f)[12]<-'f'
f<-f[12,]
f<-t(f)
save(f,file='f.RData')

#------------------------h----------------------
rownames(h)<-h$gene
comsam<-intersect(rownames(adjustedtm),rownames(h))
h<-cbind(adjustedtm[comsam,],h[,1])
write.xlsx(h,'h.xlsx')
scores <- colSums(h[, 1:ncol(h)-1] * h$coef)
h<-h[,-3047]
h<-rbind(h,scores)
rownames(h)[5]<-'h'
h<-h[5,]
h<-t(h)
save(h,file='h.RData')

#---------------p----------------------
rownames(p)<-p$gene
comsam<-intersect(rownames(adjustedtm),rownames(p))
p<-cbind(adjustedtm[comsam,],p[,1])
write.xlsx(p,'p.xlsx')
scores <- colSums(p[, 1:ncol(p)-1] * p$coef)
p<-p[,-3047]
p<-rbind(p,scores)
rownames(p)[8]<-'p'
p<-p[8,]
p<-t(p)
save(p,file='p.RData')

#----------------r------------------------
rownames(r_)<-r_$gene
comsam<-intersect(rownames(adjustedtm),rownames(r_))
r_<-cbind(adjustedtm[comsam,],r_[,1])
write.xlsx(r_,'r.xlsx')
scores <- colSums(r_[, 1:ncol(r_)-1] * r_$coef)
r_<-r_[,-3047]
r_<-rbind(r_,scores)
rownames(r_)[7]<-'r'
r_<-r_[7,]
r_<-t(r_)
save(r_,file='r.RData')

#-----------------v---------------------
rownames(v)<-v$gene
comsam<-intersect(rownames(adjustedtm),rownames(v))
v<-cbind(adjustedtm[comsam,],v[,1])
write.xlsx(v,'v.xlsx')
scores <- colSums(v[, 1:ncol(v)-1] * v$coef)
v<-v[,-3047]
v<-rbind(v,scores)
rownames(v)[4]<-'v'
v<-v[4,]
v<-t(v)
save(v,file='v.RData')

#-----------------x--------------------
rownames(x)<-x$gene
comsam<-intersect(rownames(adjustedtm),rownames(x))
x<-cbind(adjustedtm[comsam,],x[,1])
write.xlsx(x,'x.xlsx')
scores <- colSums(x[, 1:ncol(x)-1] * x$coef)
x<-x[,-3047]
x<-rbind(x,scores)
rownames(x)[11]<-'x'
x<-x[11,]
x<-t(x)
save(x,file='x.RData')

#---------------------z--------------
rownames(z)<-z$gene
comsam<-intersect(rownames(adjustedtm),rownames(z))
z<-cbind(adjustedtm[comsam,],z[,1])
write.xlsx(z,'z.xlsx')
scores <- colSums(z[, 1:ncol(z)-1] * z$coef)
z<-z[,-3047]
z<-rbind(z,scores)
rownames(z)[8]<-'z'
z<-z[8,]
z<-t(z)
save(z,file='z.RData')

#------------aa-------------------
rownames(aa)<-aa$gene
comsam<-intersect(rownames(adjustedtm),rownames(aa))
aa<-cbind(adjustedtm[comsam,],aa[,1])
write.xlsx(aa,'aa.xlsx')
scores <- colSums(aa[, 1:ncol(aa)-1] * aa$coef)
aa<-aa[,-3047]
aa<-rbind(aa,scores)
rownames(aa)[10]<-'aa'
aa<-aa[10,]
aa<-t(aa)
save(aa,file='aa.RData')

#----------------ab------------------
rownames(ab)<-ab$gene
comsam<-intersect(rownames(adjustedtm),rownames(ab))
ab<-cbind(adjustedtm[comsam,],ab[,1])
write.xlsx(ab,'ab.xlsx')
scores <- colSums(ab[, 1:ncol(ab)-1] * ab$coef)
ab<-ab[,-3047]
ab<-rbind(ab,scores)
rownames(ab)[10]<-'ab'
ab<-ab[10,]
ab<-t(ab)
save(ab,file='ab.RData')

#--------------------ac---------------
rownames(ac)<-ac$gene
comsam<-intersect(rownames(adjustedtm),rownames(ac))
ac<-cbind(adjustedtm[comsam,],ac[,1])
write.xlsx(ac,'ac.xlsx')
scores <-colSums(ac[, 1:ncol(ac)-1] * ac$coef)
ac<-as.data.frame(scores)
colnames(ac)[1]<-'ac'
save(ac,file='ac.RData')

#----------------ae------------
rownames(ae)<-ae$gene
comsam<-intersect(rownames(adjustedtm),rownames(ae))
ae<-cbind(adjustedtm[comsam,],ae[,1])
write.xlsx(ae,'ae.xlsx')
scores <- colSums(ae[, 1:ncol(ae)-1] * ae$coef)
ae<-as.data.frame(scores)
colnames(ae)[1]<-'ae'
save(ae,file='ae.RData')

#--------------af------------------
rownames(af)<-af$gene
comsam<-intersect(rownames(adjustedtm),rownames(af))
af<-cbind(adjustedtm[comsam,],af[,1])
write.xlsx(af,'af.xlsx')
scores <- colSums(af[, 1:ncol(af)-1] * af$coef)
af<-as.data.frame(scores)
colnames(af)[1]<-'af'
save(af,file='af.RData')

#----------------ah----------------
rownames(ah)<-ah$gene
comsam<-intersect(rownames(adjustedtm),rownames(ah))
ah<-cbind(adjustedtm[comsam,],ah[,1])
write.xlsx(ah,'ah.xlsx')
scores <- colSums(ah[, 1:ncol(ah)-1] * ah$coef)
ah<-as.data.frame(scores)
colnames(ah)[1]<-'ah'
save(ah,file='ah.RData')

#-----------------al---------------
rownames(al)<-al$gene
comsam<-intersect(rownames(adjustedtm),rownames(al))
al<-cbind(adjustedtm[comsam,],al[,1])
write.xlsx(al,'al.xlsx')
scores <- colSums(al[, 1:ncol(al)-1] * al$coef)
al<-as.data.frame(scores)
colnames(al)[1]<-'al'
save(al,file='al.RData')

#---------am---------------
rownames(am)<-am$gene
comsam<-intersect(rownames(adjustedtm),rownames(am))
am<-cbind(adjustedtm[comsam,],am[,1])
write.xlsx(am,'am.xlsx')
scores <- colSums(am[, 1:ncol(am)-1] * am$coef)
am<-as.data.frame(scores)
colnames(am)[1]<-'am'
save(am,file='am.RData')

#-------------aq-----------------
rownames(aq)<-aq$gene
comsam<-intersect(rownames(adjustedtm),rownames(aq))
aq<-cbind(adjustedtm[comsam,],aq[,1])
write.xlsx(aq,'aq.xlsx')
scores <- colSums(aq[, 1:ncol(aq)-1] * aq$coef)
aq<-as.data.frame(scores)
colnames(aq)[1]<-'aq'
save(aq,file='aq.RData')

#--------------as-----------------
rownames(as)<-as$gene
comsam<-intersect(rownames(adjustedtm),rownames(as))
as<-cbind(adjustedtm[comsam,],as[,1])
write.xlsx(as,'as.xlsx')
scores <- colSums(as[, 1:ncol(as)-1] * as$coef)
as<-as.data.frame(scores)
colnames(as)[1]<-'as'
save(as,file='as.RData')

#-----------aw------------
rownames(aw)<-aw$gene
comsam<-intersect(rownames(adjustedtm),rownames(aw))
aw<-cbind(adjustedtm[comsam,],aw[,1])
write.xlsx(aw,'aw.xlsx')
scores <- colSums(aw[, 1:ncol(aw)-1] * aw$coef)
aw<-as.data.frame(scores)
colnames(aw)[1]<-'aw'
save(aw,file='aw.RData')

#------ay--------
rownames(ay)<-ay$gene
comsam<-intersect(rownames(adjustedtm),rownames(ay))
ay<-cbind(adjustedtm[comsam,],ay[,1])
write.xlsx(ay,'ay.xlsx')
scores <- colSums(ay[, 1:ncol(ay)-1] * ay$coef)
ay<-as.data.frame(scores)
colnames(ay)[1]<-'ay'
save(ay,file='ay.RData')

#-------az--------------
rownames(az)<-az$gene
comsam<-intersect(rownames(adjustedtm),rownames(az))
az<-cbind(adjustedtm[comsam,],az[,1])
write.xlsx(az,'az.xlsx')
scores <- colSums(az[, 1:ncol(az)-1] * az$coef)
az<-as.data.frame(scores)
colnames(az)[1]<-'az'
save(az,file='az.RData')

#----------bg-----------
rownames(bg)<-bg$gene
comsam<-intersect(rownames(adjustedtm),rownames(bg))
bg<-cbind(adjustedtm[comsam,],bg[,1])
write.xlsx(bg,'bg.xlsx')
scores <- colSums(bg[, 1:ncol(bg)-1] * bg$coef)
bg<-as.data.frame(scores)
colnames(bg)[1]<-'bg'
save(bg,file='bg.RData')

#------bl-------------
rownames(bl)<-bl$gene
comsam<-intersect(rownames(adjustedtm),rownames(bl))
bl<-cbind(adjustedtm[comsam,],bl[,1])
write.xlsx(bl,'bl.xlsx')
scores <- colSums(bl[, 1:ncol(bl)-1] * bl$coef)
bl<-as.data.frame(scores)
colnames(bl)[1]<-'bl'
save(bl,file='bl.RData')

#------------bn-------
rownames(bn)<-bn$gene
comsam<-intersect(rownames(adjustedtm),rownames(bn))
bn<-cbind(adjustedtm[comsam,],bn[,1])
write.xlsx(bn,'bn.xlsx')
scores <- colSums(bn[, 1:ncol(bn)-1] * bn$coef)
bn<-as.data.frame(scores)
colnames(bn)[1]<-'bn'
save(bn,file='bn.RData')

#----bo---------
rownames(bo)<-bo$gene
comsam<-intersect(rownames(adjustedtm),rownames(bo))
bo<-cbind(adjustedtm[comsam,],bo[,1])
write.xlsx(bo,'bo.xlsx')
scores <- colSums(bo[, 1:ncol(bo)-1] * bo$coef)
bo<-as.data.frame(scores)
colnames(bo)[1]<-'bo'
save(bo,file='bo.RData')

#-------------bp----------
rownames(bp)<-bp$gene
comsam<-intersect(rownames(adjustedtm),rownames(bp))
bp<-cbind(adjustedtm[comsam,],bp[,1])
write.xlsx(bp,'bp.xlsx')
scores <- colSums(bp[, 1:ncol(bp)-1] * bp$coef)
bp<-as.data.frame(scores)
colnames(bp)[1]<-'bp'
save(bp,file='bp.RData')

#-----------------bq-----------
rownames(bq)<-bq$gene
comsam<-intersect(rownames(adjustedtm),rownames(bq))
bq<-cbind(adjustedtm[comsam,],bq[,1])
write.xlsx(bq,'bq.xlsx')
scores <- colSums(bq[, 1:ncol(bq)-1] * bq$coef)
bq<-as.data.frame(scores)
colnames(bq)[1]<-'bq'
save(bq,file='bq.RData')

#-----------------br-----------
rownames(br)<-br$gene
comsam<-intersect(rownames(adjustedtm),rownames(br))
br<-cbind(adjustedtm[comsam,],br[,1])
write.xlsx(br,'br.xlsx')
scores <- colSums(br[, 1:ncol(br)-1] * br$coef)
br<-as.data.frame(scores)
colnames(br)[1]<-'br'
save(br,file='br.RData')

#-----------------bs-----------
rownames(bs)<-bs$gene
comsam<-intersect(rownames(adjustedtm),rownames(bs))
bs<-cbind(adjustedtm[comsam,],bs[,1])
write.xlsx(bs,'bs.xlsx')
scores <- colSums(bs[, 1:ncol(bs)-1] * bs$coef)
bs<-as.data.frame(scores)
colnames(bs)[1]<-'bs'
save(bs,file='bs.RData')

#-----------------bt-----------
rownames(bt)<-bt$gene
comsam<-intersect(rownames(adjustedtm),rownames(bt))
bt<-cbind(adjustedtm[comsam,],bt[,1])
write.xlsx(bt,'bt.xlsx')
scores <- colSums(bt[, 1:ncol(bt)-1] * bt$coef)
bt<-as.data.frame(scores)
colnames(bt)[1]<-'bt'
save(bt,file='bt.RData')

#-----------------bu-----------
rownames(bu)<-bu$gene
comsam<-intersect(rownames(adjustedtm),rownames(bu))
bu<-cbind(adjustedtm[comsam,],bu[,1])
write.xlsx(bu,'bu.xlsx')
scores <- colSums(bu[, 1:ncol(bu)-1] * bu$coef)
bu<-as.data.frame(scores)
colnames(bu)[1]<-'bu'
save(bu,file='bu.RData')

#-----------------bv-----------
rownames(bv)<-bv$gene
comsam<-intersect(rownames(adjustedtm),rownames(bv))
bv<-cbind(adjustedtm[comsam,],bv[,1])
write.xlsx(bv,'bv.xlsx')
scores <- colSums(bv[, 1:ncol(bv)-1] * bv$coef)
bv<-as.data.frame(scores)
colnames(bv)[1]<-'bv'
save(bv,file='bv.RData')

#-----------------bw-----------
rownames(bw)<-bw$gene
comsam<-intersect(rownames(adjustedtm),rownames(bw))
bw<-cbind(adjustedtm[comsam,],bw[,1])
write.xlsx(bw,'bw.xlsx')
scores <- colSums(bw[, 1:ncol(bw)-1] * bw$coef)
bw<-as.data.frame(scores)
colnames(bw)[1]<-'bw'
save(bw,file='bw.RData')

#-----------------bx-----------
rownames(bx)<-bx$gene
comsam<-intersect(rownames(adjustedtm),rownames(bx))
bx<-cbind(adjustedtm[comsam,],bx[,1])
write.xlsx(bx,'bx.xlsx')
scores <- colSums(bx[, 1:ncol(bx)-1] * bx$coef)
bx<-as.data.frame(scores)
colnames(bx)[1]<-'bx'
save(bx,file='bx.RData')

#-----------------by-----------
rownames(by)<-by$gene
comsam<-intersect(rownames(adjustedtm),rownames(by))
by<-cbind(adjustedtm[comsam,],by[,1])
write.xlsx(by,'by.xlsx')
scores <- colSums(by[, 1:ncol(by)-1] * by$coef)
by<-as.data.frame(scores)
colnames(by)[1]<-'by'
save(by,file='by.RData')

#-----------------bz-----------
rownames(bz)<-bz$gene
comsam<-intersect(rownames(adjustedtm),rownames(bz))
bz<-cbind(adjustedtm[comsam,],bz[,1])
write.xlsx(bz,'bz.xlsx')
scores <- colSums(bz[, 1:ncol(bz)-1] * bz$coef)
bz<-as.data.frame(scores)
colnames(bz)[1]<-'bz'
save(bz,file='bz.RData')

#-----------------ca-----------
rownames(ca)<-ca$gene
comsam<-intersect(rownames(adjustedtm),rownames(ca))
ca<-cbind(adjustedtm[comsam,],ca[,1])
write.xlsx(ca,'ca.xlsx')
scores <- colSums(ca[, 1:ncol(ca)-1] * ca$coef)
ca<-as.data.frame(scores)
colnames(ca)[1]<-'ca'
save(ca,file='ca.RData')

#-----------------cb-----------
rownames(cb)<-cb$gene
comsam<-intersect(rownames(adjustedtm),rownames(cb))
cb<-cbind(adjustedtm[comsam,],cb[,1])
write.xlsx(cb,'cb.xlsx')
scores <- colSums(cb[, 1:ncol(cb)-1] * cb$coef)
cb<-as.data.frame(scores)
colnames(cb)[1]<-'cb'
save(cb,file='cb.RData')

#-----------------cc-----------
rownames(cc)<-cc$gene
comsam<-intersect(rownames(adjustedtm),rownames(cc))
cc<-cbind(adjustedtm[comsam,],cc[,1])
write.xlsx(cc,'cc.xlsx')
scores <- colSums(cc[, 1:ncol(cc)-1] * cc$coef)
cc<-as.data.frame(scores)
colnames(cc)[1]<-'cc'
save(cc,file='cc.RData')

#-----------------cd-----------
rownames(cd)<-cd$gene
comsam<-intersect(rownames(adjustedtm),rownames(cd))
cd<-cbind(adjustedtm[comsam,],cd[,1])
write.xlsx(cd,'cd.xlsx')
scores <- colSums(cd[, 1:ncol(cd)-1] * cd$coef)
cd<-as.data.frame(scores)
colnames(cd)[1]<-'cd'
save(cd,file='cd.RData')

#-----------------ce-----------
rownames(ce)<-ce$gene
comsam<-intersect(rownames(adjustedtm),rownames(ce))
ce<-cbind(adjustedtm[comsam,],ce[,1])
write.xlsx(ce,'ce.xlsx')
scores <- colSums(ce[, 1:ncol(ce)-1] * ce$coef)
ce<-as.data.frame(scores)
colnames(ce)[1]<-'ce'
save(ce,file='ce.RData')

#-----------------cf-----------
rownames(cf)<-cf$gene
comsam<-intersect(rownames(adjustedtm),rownames(cf))
cf<-cbind(adjustedtm[comsam,],cf[,1])
write.xlsx(cf,'cf.xlsx')
scores <- colSums(cf[, 1:ncol(cf)-1] * cf$coef)
cf<-as.data.frame(scores)
colnames(cf)[1]<-'cf'
save(cf,file='cf.RData')

#-----------------cg-----------
rownames(cg)<-cg$gene
comsam<-intersect(rownames(adjustedtm),rownames(cg))
cg<-cbind(adjustedtm[comsam,],cg[,1])
write.xlsx(cg,'cg.xlsx')
scores <- colSums(cg[, 1:ncol(cg)-1] * cg$coef)
cg<-as.data.frame(scores)
colnames(cg)[1]<-'cg'
save(cg,file='cg.RData')

#-----------------ch-----------
rownames(ch)<-ch$gene
comsam<-intersect(rownames(adjustedtm),rownames(ch))
ch<-cbind(adjustedtm[comsam,],ch[,1])
write.xlsx(ch,'ch.xlsx')
scores <- colSums(ch[, 1:ncol(ch)-1] * ch$coef)
ch<-as.data.frame(scores)
colnames(ch)[1]<-'ch'
save(ch,file='ch.RData')

#--------------------训练集合并-----------------
comparetrain<-cbind(score_train,ae, af, al, ah, ac, aa, ab, aq, am, bl, 
                    az, ay, b, bg, as, aw, bn, bo, bp, bq, br, 
                    bs, by, bz, bx, bw, bu, bt, bv, ce, cb, ca, 
                    cd, cc, d, e, cg, cf, ch, f, h, r_, p, v, x, z)
#----------批量单因素-------------------
library(survival)
covariates<-colnames(comparetrain[,3:ncol(comparetrain)])
trainresults <- data.frame(covariate = character(0),
                      beta = numeric(0),
                      HR=numeric(0),
                      HR.confint.lower=numeric(0),
                      HR.confint.upper=numeric(0),
                      wald.test = numeric(0),
                      p.value = numeric(0))
for (covariate in covariates) {
  formula <- as.formula(paste("Surv(OS.time, OS) ~", covariate))
  model <- coxph(formula, data = comparetrain,iter.max = 1000)
  summary <- summary(model)
  p.value <- signif(summary$wald["pvalue"], digits = 2)
  wald.test <- signif(summary$wald["test"], digits = 2)
  beta <- signif(summary$coef[1], digits = 2)
  HR<-signif(summary$conf.int[,1],3)
  HR.confint.lower <- signif(summary$conf.int[,"lower .95"], 3)
  HR.confint.upper <- signif(summary$conf.int[,"upper .95"], 3)
  trainresults <- rbind(trainresults, data.frame(covariate = covariate,
                                       beta = beta,
                                       HR=HR,
                                       HR.confint.lower=HR.confint.lower,
                                       HR.confint.upper=HR.confint.upper,
                                       wald.test= wald.test,
                                       p.value = p.value))
}
save(trainresults,file = 'trainresults.RData')
##------------tcga-------------
comparetcga<-comparetrain[c(1:1068),]
tcgaresults <- data.frame(covariate = character(0),
                           beta = numeric(0),
                           HR=numeric(0),
                           HR.confint.lower=numeric(0),
                           HR.confint.upper=numeric(0),
                           wald.test = numeric(0),
                           p.value = numeric(0))
for (covariate in covariates) {
  formula <- as.formula(paste("Surv(OS.time, OS) ~", covariate))
  model <- coxph(formula, data = comparetcga,iter.max = 1000)
  summary <- summary(model)
  p.value <- signif(summary$wald["pvalue"], digits = 2)
  wald.test <- signif(summary$wald["test"], digits = 2)
  beta <- signif(summary$coef[1], digits = 2)
  HR<-signif(summary$conf.int[,1],3)
  HR.confint.lower <- signif(summary$conf.int[,"lower .95"], 3)
  HR.confint.upper <- signif(summary$conf.int[,"upper .95"], 3)
  tcgaresults <- rbind(tcgaresults, data.frame(covariate = covariate,
                                                 beta = beta,
                                                 HR=HR,
                                                 HR.confint.lower=HR.confint.lower,
                                                 HR.confint.upper=HR.confint.upper,
                                                 wald.test= wald.test,
                                                 p.value = p.value))
}
save(tcgaresults,file = 'tcgaresults.RData')

##-----------------meta-----------------
comparemeta<-comparetrain[c(1069:3046),]
metaresults <- data.frame(covariate = character(0),
                           beta = numeric(0),
                           HR=numeric(0),
                           HR.confint.lower=numeric(0),
                           HR.confint.upper=numeric(0),
                           wald.test = numeric(0),
                           p.value = numeric(0))
for (covariate in covariates) {
  formula <- as.formula(paste("Surv(OS.time, OS) ~", covariate))
  model <- coxph(formula, data = comparemeta,iter.max = 1000)
  summary <- summary(model)
  p.value <- signif(summary$wald["pvalue"], digits = 2)
  wald.test <- signif(summary$wald["test"], digits = 2)
  beta <- signif(summary$coef[1], digits = 2)
  HR<-signif(summary$conf.int[,1],3)
  HR.confint.lower <- signif(summary$conf.int[,"lower .95"], 3)
  HR.confint.upper <- signif(summary$conf.int[,"upper .95"], 3)
  metaresults <- rbind(metaresults, data.frame(covariate = covariate,
                                               beta = beta,
                                               HR=HR,
                                               HR.confint.lower=HR.confint.lower,
                                               HR.confint.upper=HR.confint.upper,
                                               wald.test= wald.test,
                                               p.value = p.value))
}
save(metaresults,file = 'metaresults.RData')

#-----------------------------GSE42568-----------------------------------------
comsam<-intersect(colnames(exprset_unique),rownames(testos))
exprset_unique<-exprset_unique[,comsam]
##---aa---
rownames(aa)<-aa$gene
comsam<-intersect(rownames(exprset_unique),rownames(aa))
aa<-cbind(exprset_unique[comsam,],aa[,1])
scores <- colSums(aa[, 1:ncol(aa)-1] * aa$coef)
aa42568<-as.data.frame(scores)
colnames(aa42568)[1]<-'aa'
save(aa42568,file='aa42568.RData')
##---ab---
rownames(ab)<-ab$gene
comsam<-intersect(rownames(exprset_unique),rownames(ab))
ab<-cbind(exprset_unique[comsam,],ab[,1])
scores <- colSums(ab[, 1:ncol(ab)-1] * ab$coef)
ab42568<-as.data.frame(scores)
colnames(ab42568)[1]<-'ab'
save(ab42568,file='ab42568.RData')
##----ac-----
rownames(ac)<-ac$gene
comsam<-intersect(rownames(exprset_unique),rownames(ac))
ac<-cbind(exprset_unique[comsam,],ac[,1])
scores <- colSums(ac[, 1:ncol(ac)-1] * ac$coef)
ac42568<-as.data.frame(scores)
colnames(ac42568)[1]<-'ac'
save(ac42568,file='ac42568.RData')
##---ae---
rownames(ae)<-ae$gene
comsam<-intersect(rownames(exprset_unique),rownames(ae))
ae<-cbind(exprset_unique[comsam,],ae[,1])
scores <- colSums(ae[, 1:ncol(ae)-1] * ae$coef)
ae42568<-as.data.frame(scores)
colnames(ae42568)[1]<-'ae'
save(ae42568,file='ae42568.RData')
##---af---
rownames(af)<-af$gene
comsam<-intersect(rownames(exprset_unique),rownames(af))
af<-cbind(exprset_unique[comsam,],af[,1])
scores <- colSums(af[, 1:ncol(af)-1] * af$coef)
af42568<-as.data.frame(scores)
colnames(af42568)[1]<-'af'
save(af42568,file='af42568.RData')
##-----ah---------
rownames(ah)<-ah$gene
comsam<-intersect(rownames(exprset_unique),rownames(ah))
ah<-cbind(exprset_unique[comsam,],ah[,1])
scores <- colSums(ah[, 1:ncol(ah)-1] * ah$coef)
ah42568<-as.data.frame(scores)
colnames(ah42568)[1]<-'ah'
save(ah42568,file='ah42568.RData')
##-----aj------
rownames(aj)<-aj$gene
comsam<-intersect(rownames(exprset_unique),rownames(aj))
aj<-cbind(exprset_unique[comsam,],aj[,1])
scores <- colSums(aj[, 1:ncol(aj)-1] * aj$coef)
aj42568<-as.data.frame(scores)
colnames(aj42568)[1]<-'aj'
save(aj42568,file='aj42568.RData')
##----al-----
rownames(al)<-al$gene
comsam<-intersect(rownames(exprset_unique),rownames(al))
al<-cbind(exprset_unique[comsam,],al[,1])
scores <- colSums(al[, 1:ncol(al)-1] * al$coef)
al42568<-as.data.frame(scores)
colnames(al42568)[1]<-'al'
save(al42568,file='al42568.RData')
##----am-----
rownames(am)<-am$gene
comsam<-intersect(rownames(exprset_unique),rownames(am))
am<-cbind(exprset_unique[comsam,],am[,1])
scores <- colSums(am[, 1:ncol(am)-1] * am$coef)
am42568<-as.data.frame(scores)
colnames(am42568)[1]<-'am'
save(am42568,file='am42568.RData')
##------aq-----
rownames(aq)<-aq$gene
comsam<-intersect(rownames(exprset_unique),rownames(aq))
aq<-cbind(exprset_unique[comsam,],aq[,1])
scores <- colSums(aq[, 1:ncol(aq)-1] * aq$coef)
aq42568<-as.data.frame(scores)
colnames(aq42568)[1]<-'aq'
save(aq42568,file='aq42568.RData')
##----ar-----
rownames(ar)<-ar$gene
comsam<-intersect(rownames(exprset_unique),rownames(ar))
ar<-cbind(exprset_unique[comsam,],ar[,1])
scores <- colSums(ar[, 1:ncol(ar)-1] * ar$coef)
ar42568<-as.data.frame(scores)
colnames(ar42568)[1]<-'ar'
save(ar42568,file='ar42568.RData')
##----as-----
rownames(as)<-as$gene
comsam<-intersect(rownames(exprset_unique),rownames(as))
as<-cbind(exprset_unique[comsam,],as[,1])
scores <- colSums(as[, 1:ncol(as)-1] * as$coef)
as42568<-as.data.frame(scores)
colnames(as42568)[1]<-'as'
save(as42568,file='as42568.RData')
##----aw-----
rownames(aw)<-aw$gene
comsam<-intersect(rownames(exprset_unique),rownames(aw))
aw<-cbind(exprset_unique[comsam,],aw[,1])
scores <- colSums(aw[, 1:ncol(aw)-1] * aw$coef)
aw42568<-as.data.frame(scores)
colnames(aw42568)[1]<-'aw'
save(aw42568,file='aw42568.RData')
##----ay-----
rownames(ay)<-ay$gene
comsam<-intersect(rownames(exprset_unique),rownames(ay))
ay<-cbind(exprset_unique[comsam,],ay[,1])
scores <- colSums(ay[, 1:ncol(ay)-1] * ay$coef)
ay42568<-as.data.frame(scores)
colnames(ay42568)[1]<-'ay'
save(ay42568,file='ay42568.RData')
##-----az-----
rownames(az)<-az$gene
comsam<-intersect(rownames(exprset_unique),rownames(az))
az<-cbind(exprset_unique[comsam,],az[,1])
scores <- colSums(az[, 1:ncol(az)-1] * az$coef)
az42568<-as.data.frame(scores)
colnames(az42568)[1]<-'az'
save(az42568,file='az42568.RData')
##------b-------------
rownames(b)<-b$gene
comsam<-intersect(rownames(exprset_unique),rownames(b))
b<-cbind(exprset_unique[comsam,],b[,1])
scores <- colSums(b[, 1:ncol(b)-1] * b$coef)
b42568<-as.data.frame(scores)
colnames(b42568)[1]<-'b'
save(b42568,file='b42568.RData')
##----bg-----
rownames(bg)<-bg$gene
comsam<-intersect(rownames(exprset_unique),rownames(bg))
bg<-cbind(exprset_unique[comsam,],bg[,1])
scores <- colSums(bg[, 1:ncol(bg)-1] * bg$coef)
bg42568<-as.data.frame(scores)
colnames(bg42568)[1]<-'bg'
save(bg42568,file='bg42568.RData')
##----bn-----
rownames(bn)<-bn$gene
comsam<-intersect(rownames(exprset_unique),rownames(bn))
bn<-cbind(exprset_unique[comsam,],bn[,1])
scores <- colSums(bn[, 1:ncol(bn)-1] * bn$coef)
bn42568<-as.data.frame(scores)
colnames(bn42568)[1]<-'bn'
save(bn42568,file='bn42568.RData')
##-----bo-----
rownames(bo)<-bo$gene
comsam<-intersect(rownames(exprset_unique),rownames(bo))
bo<-cbind(exprset_unique[comsam,],bo[,1])
scores <- colSums(bo[, 1:ncol(bo)-1] * bo$coef)
bo42568<-as.data.frame(scores)
colnames(bo42568)[1]<-'bo'
save(bo42568,file='bo42568.RData')
##----bp-----
rownames(bp)<-bp$gene
comsam<-intersect(rownames(exprset_unique),rownames(bp))
bp<-cbind(exprset_unique[comsam,],bp[,1])
scores <- colSums(bp[, 1:ncol(bp)-1] * bp$coef)
bp20685<-as.data.frame(scores)
colnames(bp20685)[1]<-'bp'
save(bp20685,file='bp20685.RData')
##----bq-----
rownames(bq)<-bq$gene
comsam<-intersect(rownames(exprset_unique),rownames(bq))
bq<-cbind(exprset_unique[comsam,],bq[,1])
scores <- colSums(bq[, 1:ncol(bq)-1] * bq$coef)
bq42568<-as.data.frame(scores)
colnames(bq42568)[1]<-'bq'
save(bq42568,file='bq42568.RData')
##---br-----
rownames(br)<-br$gene
comsam<-intersect(rownames(exprset_unique),rownames(br))
br<-cbind(exprset_unique[comsam,],br[,1])
scores <- colSums(br[, 1:ncol(br)-1] * br$coef)
br42568<-as.data.frame(scores)
colnames(br42568)[1]<-'br'
save(br42568,file='br42568.RData')
##----bs----
rownames(bs)<-bs$gene
comsam<-intersect(rownames(exprset_unique),rownames(bs))
bs<-cbind(exprset_unique[comsam,],bs[,1])
scores <- colSums(bs[, 1:ncol(bs)-1] * bs$coef)
bs42568<-as.data.frame(scores)
colnames(bs42568)[1]<-'bs'
save(bs42568,file='bs42568.RData')
##----bt-----
#comsam不匹配

##----bu-----
rownames(bu)<-bu$gene
comsam<-intersect(rownames(exprset_unique),rownames(bu))
bu<-cbind(exprset_unique[comsam,],bu[,1])
scores <- colSums(bu[, 1:ncol(bu)-1] * bu$coef)
bu42568<-as.data.frame(scores)
colnames(bu42568)[1]<-'bu'
save(bu42568,file='bu42568.RData')
##----bv----
rownames(bv)<-bv$gene
comsam<-intersect(rownames(exprset_unique),rownames(bv))
bv<-cbind(exprset_unique[comsam,],bv[,1])
scores <- colSums(bv[, 1:ncol(bv)-1] * bv$coef)
bv42568<-as.data.frame(scores)
colnames(bv42568)[1]<-'bv'
save(bv42568,file='bv42568.RData')
##-----bw----
rownames(bw)<-bw$gene
comsam<-intersect(rownames(exprset_unique),rownames(bw))
bw<-cbind(exprset_unique[comsam,],bw[,1])
scores <- colSums(bw[, 1:ncol(bw)-1] * bw$coef)
bw42568<-as.data.frame(scores)
colnames(bw42568)[1]<-'bw'
save(bw42568,file='bw42568.RData')
##----bx-----
rownames(bx)<-bx$gene
comsam<-intersect(rownames(exprset_unique),rownames(bx))
bx<-cbind(exprset_unique[comsam,],bx[,1])
scores <- colSums(bx[, 1:ncol(bx)-1] * bx$coef)
bx42568<-as.data.frame(scores)
colnames(bx42568)[1]<-'bx'
save(bx42568,file='bx42568.RData')
##----by-----
rownames(by)<-by$gene
comsam<-intersect(rownames(exprset_unique),rownames(by))
by<-cbind(exprset_unique[comsam,],by[,1])
scores <- colSums(by[, 1:ncol(by)-1] * by$coef)
by42568<-as.data.frame(scores)
colnames(by42568)[1]<-'by'
save(by42568,file='by42568.RData')
##----bz-----
rownames(bz)<-bz$gene
comsam<-intersect(rownames(exprset_unique),rownames(bz))
bz<-cbind(exprset_unique[comsam,],bz[,1])
scores <- colSums(bz[, 1:ncol(bz)-1] * bz$coef)
bz42568<-as.data.frame(scores)
colnames(bz42568)[1]<-'bz'
save(bz42568,file='bz42568.RData')
##----ca-----
#comsam不同
##----cb-----
rownames(cb)<-cb$gene
comsam<-intersect(rownames(exprset_unique),rownames(cb))
cb<-cbind(exprset_unique[comsam,],cb[,1])
scores <- colSums(cb[, 1:ncol(cb)-1] * cb$coef)
cb42568<-as.data.frame(scores)
colnames(cb42568)[1]<-'cb'
save(cb42568,file='cb42568.RData')
##----cc-----
rownames(cc)<-cc$gene
comsam<-intersect(rownames(exprset_unique),rownames(cc))
cc<-cbind(exprset_unique[comsam,],cc[,1])
scores <- colSums(cc[, 1:ncol(cc)-1] * cc$coef)
cc42568<-as.data.frame(scores)
colnames(cc42568)[1]<-'cc'
save(cc42568,file='cc42568.RData')
##----cd-----
#comsam不同
##----ce-----
rownames(ce)<-ce$gene
comsam<-intersect(rownames(exprset_unique),rownames(ce))
ce<-cbind(exprset_unique[comsam,],ce[,1])
scores <- colSums(ce[, 1:ncol(ce)-1] * ce$coef)
ce42568<-as.data.frame(scores)
colnames(ce42568)[1]<-'ce'
save(ce42568,file='ce42568.RData')
##----cf-----
rownames(cf)<-cf$gene
comsam<-intersect(rownames(exprset_unique),rownames(cf))
cf<-cbind(exprset_unique[comsam,],cf[,1])
scores <- colSums(cf[, 1:ncol(cf)-1] * cf$coef)
cf42568<-as.data.frame(scores)
colnames(cf42568)[1]<-'cf'
save(cf42568,file='cf42568.RData')
##----cg-----
rownames(cg)<-cg$gene
comsam<-intersect(rownames(exprset_unique),rownames(cg))
cg<-cbind(exprset_unique[comsam,],cg[,1])
scores <- colSums(cg[, 1:ncol(cg)-1] * cg$coef)
cg42568<-as.data.frame(scores)
colnames(cg42568)[1]<-'cg'
save(cg42568,file='cg42568.RData')
##----ch-----
rownames(ch)<-ch$gene
comsam<-intersect(rownames(exprset_unique),rownames(ch))
ch<-cbind(exprset_unique[comsam,],ch[,1])
scores <- colSums(ch[, 1:ncol(ch)-1] * ch$coef)
ch42568<-as.data.frame(scores)
colnames(ch42568)[1]<-'ch'
save(ch42568,file='ch42568.RData')
##-----d-----
rownames(d)<-d$gene
comsam<-intersect(rownames(exprset_unique),rownames(d))
d<-cbind(exprset_unique[comsam,],d[,1])
scores <- colSums(d[, 1:ncol(d)-1] * d$coef)
d42568<-as.data.frame(scores)
colnames(d42568)[1]<-'d'
save(d42568,file='d42568.RData')
##----e-----
rownames(e)<-e$gene
comsam<-intersect(rownames(exprset_unique),rownames(e))
e<-cbind(exprset_unique[comsam,],e[,1])
scores <- colSums(e[, 1:ncol(e)-1] * e$coef)
e42568<-as.data.frame(scores)
colnames(e42568)[1]<-'e'
save(e42568,file='e42568.RData')
##-----f-----
rownames(f)<-f$gene
comsam<-intersect(rownames(exprset_unique),rownames(f))
f<-cbind(exprset_unique[comsam,],f[,1])
scores <- colSums(f[, 1:ncol(f)-1] * f$coef)
f42568<-as.data.frame(scores)
colnames(f42568)[1]<-'f'
save(f42568,file='f42568.RData')
##----h-----
rownames(h)<-h$gene
comsam<-intersect(rownames(exprset_unique),rownames(h))
h<-cbind(exprset_unique[comsam,],h[,1])
scores <- colSums(h[, 1:ncol(h)-1] * h$coef)
h42568<-as.data.frame(scores)
colnames(h42568)[1]<-'h'
save(h42568,file='h42568.RData')
##----p----
rownames(p)<-p$gene
comsam<-intersect(rownames(exprset_unique),rownames(p))
p<-cbind(exprset_unique[comsam,],p[,1])
scores <- colSums(p[, 1:ncol(p)-1] * p$coef)
p42568<-as.data.frame(scores)
colnames(p42568)[1]<-'p'
save(p42568,file='p42568.RData')
##----r----
rownames(r_)<-r_$gene
comsam<-intersect(rownames(exprset_unique),rownames(r_))
r<-cbind(exprset_unique[comsam,],r_[,1])
scores <- colSums(r[, 1:ncol(r)-1] * r$coef)
r42568<-as.data.frame(scores)
colnames(r42568)[1]<-'r'
save(r42568,file='r42568.RData')
##-----v----
rownames(v)<-v$gene
comsam<-intersect(rownames(exprset_unique),rownames(v))
v<-cbind(exprset_unique[comsam,],v[,1])
scores <- colSums(v[, 1:ncol(v)-1] * v$coef)
v42568<-as.data.frame(scores)
colnames(v42568)[1]<-'v'
save(v42568,file='v42568.RData')
##-----x----
rownames(x)<-x$gene
comsam<-intersect(rownames(exprset_unique),rownames(x))
x<-cbind(exprset_unique[comsam,],x[,1])
scores <- colSums(x[, 1:ncol(x)-1] * x$coef)
x42568<-as.data.frame(scores)
colnames(x42568)[1]<-'x'
save(x42568,file='x42568.RData')
##------z-----
rownames(z)<-z$gene
comsam<-intersect(rownames(exprset_unique),rownames(z))
z<-cbind(exprset_unique[comsam,],z[,1])
scores <- colSums(z[, 1:ncol(z)-1] * z$coef)
z42568<-as.data.frame(scores)
colnames(z42568)[1]<-'z'
save(z42568,file='z42568.RData')

#-----------------------42568-------------
comsam<-intersect(rownames(score_test),rownames(as42568))
score_42568<-score_test[comsam,]
compare42568<- cbind(score_42568,as42568, am42568, aj42568, aq42568, ar42568, al42568, 
                       ah42568, af42568,  aa42568, ab42568, 
                       aw42568,  ac42568, bp42568, bo42568, br42568, 
                       bq42568, bs42568,bn42568, bg42568, bx42568, 
                       b42568, bu42568, bw42568, bz42568, cc42568, cb42568, 
                       h42568, by42568, e42568, f42568, cg42568, d42568, 
                       ch42568, cf42568, ce42568, x42568, v42568, 
                        z42568, p42568, r42568 )
save(compare42568,file = 'compare42568.RData')
library(survival)
covariates<-colnames(compare42568[,3:ncol(compare42568)])
results42568 <- data.frame(covariate = character(0),
                           beta = numeric(0),
                           HR=numeric(0),
                           HR.confint.lower=numeric(0),
                           HR.confint.upper=numeric(0),
                           wald.test = numeric(0),
                           p.value = numeric(0))
for (covariate in covariates) {
  formula <- as.formula(paste("Surv(OS.time, OS) ~", covariate))
  model <- coxph(formula, data = compare42568,iter.max = 1000)
  summary <- summary(model)
  p.value <- signif(summary$wald["pvalue"], digits = 2)
  wald.test <- signif(summary$wald["test"], digits = 2)
  beta <- signif(summary$coef[1], digits = 2)
  HR<-signif(summary$conf.int[,1],3)
  HR.confint.lower <- signif(summary$conf.int[,"lower .95"], 3)
  HR.confint.upper <- signif(summary$conf.int[,"upper .95"], 3)
  results42568 <- rbind(results42568, data.frame(covariate = covariate,
                                                 beta = beta,
                                                 HR=HR,
                                                 HR.confint.lower=HR.confint.lower,
                                                 HR.confint.upper=HR.confint.upper,
                                                 wald.test= wald.test,
                                                 p.value = p.value))
}
save(results42568,file='results42568.RData')

#--------88770------------
comsam<-intersect(rownames(score_test),rownames(aa88770))
score_88770<-score_test[comsam,]
compare88770<- cbind(score_88770,ac88770, ab88770, 
                       aq88770, ar88770, aa88770, am88770, al88770, aj88770, 
                       ah88770, af88770, bo88770, bp88770, as88770, aw88770, 
                       b88770,  ce88770, by88770, bx88770, cc88770, cb88770, 
                       bz88770,  bw88770, br88770, bu88770, bs88770, bq88770, 
                       h88770, r88770, v88770, p88770, f88770, d88770, e88770, 
                       ch88770, cg88770,  cf88770, z88770, x88770,bg88770,bn88770)
save(compare88770,file = 'compare88770.RData')
covariates<-colnames(compare88770[,3:ncol(compare88770)])
results88770 <- data.frame(covariate = character(0),
                           beta = numeric(0),
                           HR=numeric(0),
                           HR.confint.lower=numeric(0),
                           HR.confint.upper=numeric(0),
                           wald.test = numeric(0),
                           p.value = numeric(0))
for (covariate in covariates) {
  formula <- as.formula(paste("Surv(OS.time, OS) ~", covariate))
  model <- coxph(formula, data = compare88770,iter.max = 1000)
  summary <- summary(model)
  p.value <- signif(summary$wald["pvalue"], digits = 2)
  wald.test <- signif(summary$wald["test"], digits = 2)
  beta <- signif(summary$coef[1], digits = 2)
  HR<-signif(summary$conf.int[,1],3)
  HR.confint.lower <- signif(summary$conf.int[,"lower .95"], 3)
  HR.confint.upper <- signif(summary$conf.int[,"upper .95"], 3)
  results88770 <- rbind(results88770, data.frame(covariate = covariate,
                                                 beta = beta,
                                                 HR=HR,
                                                 HR.confint.lower=HR.confint.lower,
                                                 HR.confint.upper=HR.confint.upper,
                                                 wald.test= wald.test,
                                                 p.value = p.value))
}
save(results88770,file='results88770.RData')

#-------------48390--------------
comsam<-intersect(rownames(score_test),rownames(aa48390))
score_48390<-score_test[comsam,]
compare48390<- cbind(score_48390,aw48390, ah48390, aj48390, as48390, am48390, 
                       al48390, ar48390, aq48390, af48390, ae48390, 
                       ac48390, aa48390, ay48390, ab48390, bx48390, 
                       bq48390, bw48390, bp48390, bv48390, bu48390, 
                       bs48390, br48390, bo48390, bn48390, b48390, 
                       bg48390, az48390, cb48390, ce48390, cf48390, 
                       cc48390, cg48390, by48390, bz48390, h48390, 
                       d48390, f48390, ch48390, e48390, r48390, 
                       p48390, v48390, z48390, x48390)
covariates<-colnames(compare48390[,3:ncol(compare48390)])
results48390 <- data.frame(covariate = character(0),
                           beta = numeric(0),
                           HR=numeric(0),
                           HR.confint.lower=numeric(0),
                           HR.confint.upper=numeric(0),
                           wald.test = numeric(0),
                           p.value = numeric(0))
for (covariate in covariates) {
  formula <- as.formula(paste("Surv(OS.time, OS) ~", covariate))
  model <- coxph(formula, data = compare48390,iter.max = 1000)
  summary <- summary(model)
  p.value <- signif(summary$wald["pvalue"], digits = 2)
  wald.test <- signif(summary$wald["test"], digits = 2)
  beta <- signif(summary$coef[1], digits = 2)
  HR<-signif(summary$conf.int[,1],3)
  HR.confint.lower <- signif(summary$conf.int[,"lower .95"], 3)
  HR.confint.upper <- signif(summary$conf.int[,"upper .95"], 3)
  results48390 <- rbind(results48390, data.frame(covariate = covariate,
                                                 beta = beta,
                                                 HR=HR,
                                                 HR.confint.lower=HR.confint.lower,
                                                 HR.confint.upper=HR.confint.upper,
                                                 wald.test= wald.test,
                                                 p.value = p.value))
}
save(results48390,file='results48390.RData')

#-----------------20685---------------------
comsam<-intersect(rownames(score_test),rownames(aa20685))
score_20685<-score_test[comsam,]
compare20685<- cbind(score_20685,v20685, e20685, r20685, z20685, f20685,
                       h20685, p20685, x20685, cg20685, ch20685, 
                       ce20685, cf20685, d20685, bz20685, cc20685, 
                       bx20685, bq20685, by20685, bn20685, cb20685, 
                       bo20685, bs20685, bu20685, bw20685, 
                       br20685, aq20685, am20685, ar20685,  
                       as20685, bg20685, aw20685,  b20685, 
                       ah20685, aj20685, af20685, al20685, 
                       aa20685, ac20685, ab20685,bp20685)
save(compare20685,file = 'compare20685.RData')
covariates<-colnames(compare20685[,3:ncol(compare20685)])
results20685 <- data.frame(covariate = character(0),
                           beta = numeric(0),
                           HR=numeric(0),
                           HR.confint.lower=numeric(0),
                           HR.confint.upper=numeric(0),
                           wald.test = numeric(0),
                           p.value = numeric(0))
for (covariate in covariates) {
  formula <- as.formula(paste("Surv(OS.time, OS) ~", covariate))
  model <- coxph(formula, data = compare20685,iter.max = 1000)
  summary <- summary(model)
  p.value <- signif(summary$wald["pvalue"], digits = 2)
  wald.test <- signif(summary$wald["test"], digits = 2)
  beta <- signif(summary$coef[1], digits = 2)
  HR<-signif(summary$conf.int[,1],3)
  HR.confint.lower <- signif(summary$conf.int[,"lower .95"], 3)
  HR.confint.upper <- signif(summary$conf.int[,"upper .95"], 3)
  results20685 <- rbind(results20685, data.frame(covariate = covariate,
                                                 beta = beta,
                                                 HR=HR,
                                                 HR.confint.lower=HR.confint.lower,
                                                 HR.confint.upper=HR.confint.upper,
                                                 wald.test= wald.test,
                                                 p.value = p.value))
}
save(results20685,file='results20685.RData')

#------合并-----
rownames(trainresults)<-trainresults$covariate
rownames(results88770)<-results88770$covariate
rownames(results42568)<-results42568$covariate
rownames(results20685)<-results20685$covariate
rownames(tcgaresults)<-tcgaresults$covariate
rownames(metaresults)<-metaresults$covariate
rownames(results48390)<-results48390$covariate
comsam<-intersect(rownames(trainresults),rownames(results88770))
trainresults<-trainresults[comsam,]
results20685<-results20685[comsam,]
results42568<-results42568[comsam,]
results88770<-results88770[comsam,]
results48390<-results48390[comsam,]
tcgaresults<-tcgaresults[comsam,]
metaresults<-metaresults[comsam,]
colnames(trainresults)[7]<-'train.p'
colnames(results20685)[7]<-'20685.p'
colnames(results42568)[7]<-'42568.p'
colnames(results88770)[7]<-'88770.p'
colnames(tcgaresults)[7]<-'tcga.p'
colnames(metaresults)[7]<-'meta.p'
colnames(results48390)[7]<-'48390.p'
train.p <- trainresults[, 7]
p20685 <- results20685[, 7]
p42568 <- results42568[, 7]
p88770 <- results88770[, 7]
ptcga <- tcgaresults[, 7]
pmeta <- metaresults[, 7]
p48390<-results48390[,7]
hrtrain<-trainresults[,3]
hr20685<-results20685[,3]
hr42568<-results42568[,3]
hr88770<-results88770[,3]
hrtcga<-tcgaresults[,3]
hrmeta<-metaresults[,3]
hr48390<-results48390[,3]
heatp<-cbind(train.p,ptcga,pmeta,p20685,
             p42568,p48390,p88770)
hr<-cbind(hrtrain,hrtcga,hrmeta,hr20685,
          hr42568,hr48390,hr88770)
rownames(heatp)<-rownames(trainresults)
heatp<-as.data.frame(heatp)
rownames(hr)<-rownames(trainresults)
hr<-as.data.frame(hr)
write.xlsx(heatp,'heatp2.xlsx')
#----plot-----
library(ggplot2)
library(reshape2)
hr$Variable<-rownames(hr)
# 使用melt函数融化数据框，指定Variable列作为id变量
df_long <- melt(heatp, id.vars = "Variable", variable.name = "DataSet", value.name = "PValue")
hr_df_long <- melt(hr, id.vars = "Variable",variable.name = "DataSet", value.name = "HR")
df_long <- df_long[order(df_long$Variable == 'Score', decreasing = TRUE),]
hr_df_long <- hr_df_long[order(hr_df_long$Variable == 'Score', decreasing = TRUE),]
# 将HR值合并到P值的数据框中
df_long$HR <- hr_df_long$HR
# 计算每个观测的颜色
df_long$Color <- ifelse(df_long$PValue > 0.05, 'p>0.05',
                        ifelse(df_long$HR > 1, "risky", "protective"))
df_long$Variable <- factor(df_long$Variable, levels = c("Score", setdiff(df_long$Variable, "Score")))
save(df_long,file = '模型比较长数据.RData')
p <- ggplot(df_long, aes(x = Variable, y = DataSet, fill = Color)) +
  geom_tile(color = "grey") +
  scale_fill_manual(values = c("risky" = "#E64B35", "protective" = "#8491B4B2", 'p>0.05' = '#FFFFFF')) +
  theme_minimal() +
  theme(
    text = element_text(size = 10, family = "Arial"),  # 调整文本大小和字体
    axis.title = element_text(size = 12, face = "bold"),  # 调整轴标题大小
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # 调整x轴标签角度和大小
    axis.text.y = element_text(size = 10),  # 调整y轴标签大小
    legend.title = element_text(size = 10),  # 调整图例标题大小
    legend.text = element_text(size = 8),  # 调整图例文本大小
    panel.grid.major = element_blank(),  # 去除主要网格线
    panel.grid.minor = element_blank(),  # 去除次要网格线
    panel.background = element_blank(),  # 去除背景
    legend.position = "right"  # 将图例位置设为右侧
  ) +
  labs(x = "变量", y = "数据集", fill = "分类") +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
  coord_fixed(ratio = 1) # 确保每个tile为正方形

# 调整图形尺寸
p = p + theme(
  legend.text=element_text(size=8),  # 调整图例文本大小
  legend.title=element_text(size=8),  # 调整图例标题大小
  legend.key.size = unit(0.2, "cm")  # 调整图例图标大小
)
# 显示绘图
print(p)

#-----------Cindex-------------
##-----train-----
set.seed(23)
library(survival)
# 创建一个空的 dataframe 存储结果
traincindex <- data.frame(variable = character(), cindex = numeric(), lower_ci = numeric(), upper_ci = numeric(), p_value = numeric())
for (i in 3:ncol(comparetrain) {
  fit <- coxph(Surv(OS.time, OS) ~ comparetrain[, i], data = comparetrain)
  sum.surv <- summary(fit)
  c_index <- sum.surv$concordance
  ci <- confint(fit)
  lower_ci <- ci[1]
  upper_ci <- ci[2]
  p_value <- sum.surv$coefficients[,"Pr(>|z|)"]
  traincindex <- rbind(traincindex, data.frame(variable = colnames(comparetcga)[i], 
                                               cindex = c_index, 
                                               lower_ci = lower_ci, 
                                               upper_ci = upper_ci,
                                               p_value = p_value))
}
traincindex <- traincindex %>%
  group_by(variable) %>%
  summarise(
    C_index_sum = sum(cindex),
    CI_lower_avg = mean(lower_ci),
    CI_upper_avg = mean(upper_ci),
    p_value_avg = mean(p_value)
  )
save(traincindex, file = 'traincindex.RData')
##-----tcga-----
tcgacindex <- data.frame(variable = character(), cindex = numeric(), lower_ci = numeric(), upper_ci = numeric(), p_value = numeric())
for (i in 3:ncol(comparetcga)) {
  fit <- coxph(Surv(OS.time, OS) ~ comparetcga[, i], data = comparetcga)
  sum.surv <- summary(fit)
  c_index <- sum.surv$concordance
  ci <- confint(fit)
  lower_ci <- ci[1]
  upper_ci <- ci[2]
  p_value <- sum.surv$coefficients[,"Pr(>|z|)"]
  tcgacindex <- rbind(tcgacindex, data.frame(variable = colnames(comparetcga)[i], 
                                               cindex = c_index, 
                                               lower_ci = lower_ci, 
                                               upper_ci = upper_ci,
                                               p_value = p_value))
}
tcgacindex <- tcgacindex %>%
  group_by(variable) %>%
  summarise(
    C_index_sum = sum(cindex),
    CI_lower_avg = mean(lower_ci),
    CI_upper_avg = mean(upper_ci),
    p_value_avg = mean(p_value)
  )
save(tcgacindex, file = 'tcgacindex.RData')
##-----meta-----
metacindex <- data.frame(variable = character(), cindex = numeric(), lower_ci = numeric(), upper_ci = numeric(), p_value = numeric())
for (i in 3:ncol(comparemeta)) {
  fit <- coxph(Surv(OS.time, OS) ~ comparemeta[, i], data = comparemeta)
  sum.surv <- summary(fit)
  c_index <- sum.surv$concordance
  ci <- confint(fit)
  lower_ci <- ci[1]
  upper_ci <- ci[2]
  p_value <- sum.surv$coefficients[,"Pr(>|z|)"]
  metacindex <- rbind(metacindex, data.frame(variable = colnames(comparemeta)[i], 
                                             cindex = c_index, 
                                             lower_ci = lower_ci, 
                                             upper_ci = upper_ci,
                                             p_value = p_value))
}
metacindex <- metacindex %>%
  group_by(variable) %>%
  summarise(
    C_index_sum = sum(cindex),
    CI_lower_avg = mean(lower_ci),
    CI_upper_avg = mean(upper_ci),
    p_value_avg = mean(p_value)
  )
save(metacindex, file = 'metacindex.RData')
##-----20685----
cindex20685 <- data.frame(variable = character(), cindex = numeric(), lower_ci = numeric(), upper_ci = numeric(), p_value = numeric())
for (i in 3:ncol(compare20685)) {
  fit <- coxph(Surv(OS.time, OS) ~ compare20685[, i], data = compare20685)
  sum.surv <- summary(fit)
  c_index <- sum.surv$concordance
  ci <- confint(fit)
  lower_ci <- ci[1]
  upper_ci <- ci[2]
  p_value <- sum.surv$coefficients[,"Pr(>|z|)"]
  cindex20685 <- rbind(cindex20685, data.frame(variable = colnames(compare20685)[i], 
                                             cindex = c_index, 
                                             lower_ci = lower_ci, 
                                             upper_ci = upper_ci,
                                             p_value = p_value))
}
cindex20685 <- cindex20685 %>%
  group_by(variable) %>%
  summarise(
    C_index_sum = sum(cindex),
    CI_lower_avg = mean(lower_ci),
    CI_upper_avg = mean(upper_ci),
    p_value_avg = mean(p_value)
  )
save(cindex20685,file = 'cindex20685.RData')
##-----42568----
cindex42568 <- data.frame(variable = character(), cindex = numeric(), lower_ci = numeric(), upper_ci = numeric(), p_value = numeric())
for (i in 3:ncol(compare42568)) {
  fit <- coxph(Surv(OS.time, OS) ~ compare42568[, i], data = compare42568)
  sum.surv <- summary(fit)
  c_index <- sum.surv$concordance
  ci <- confint(fit)
  lower_ci <- ci[1]
  upper_ci <- ci[2]
  p_value <- sum.surv$coefficients[,"Pr(>|z|)"]
  cindex42568 <- rbind(cindex42568, data.frame(variable = colnames(compare42568)[i], 
                                               cindex = c_index, 
                                               lower_ci = lower_ci, 
                                               upper_ci = upper_ci,
                                               p_value = p_value))
}
cindex42568 <- cindex42568 %>%
  group_by(variable) %>%
  summarise(
    C_index_sum = sum(cindex),
    CI_lower_avg = mean(lower_ci),
    CI_upper_avg = mean(upper_ci),
    p_value_avg = mean(p_value)
  )
save(cindex42568,file = 'cindex42568.RData')
##--------cindex48390---------
cindex48390 <- data.frame(variable = character(), cindex = numeric(), lower_ci = numeric(), upper_ci = numeric(), p_value = numeric())
for (i in 3:ncol(compare20685)) {
  fit <- coxph(Surv(OS.time, OS) ~ compare20685[, i], data = compare20685)
  sum.surv <- summary(fit)
  c_index <- sum.surv$concordance
  ci <- confint(fit)
  lower_ci <- ci[1]
  upper_ci <- ci[2]
  p_value <- sum.surv$coefficients[,"Pr(>|z|)"]
  cindex48390 <- rbind(cindex48390, data.frame(variable = colnames(compare20685)[i], 
                                               cindex = c_index, 
                                               lower_ci = lower_ci, 
                                               upper_ci = upper_ci,
                                               p_value = p_value))
}
cindex48390 <- cindex48390 %>%
  group_by(variable) %>%
  summarise(
    C_index_sum = sum(cindex),
    CI_lower_avg = mean(lower_ci),
    CI_upper_avg = mean(upper_ci),
    p_value_avg = mean(p_value)
  )
save(cindex48390,file = 'cindex48390.RData')


##----88770----
cindex88770 <- data.frame(variable = character(), cindex = numeric(), lower_ci = numeric(), upper_ci = numeric(), p_value = numeric())
for (i in 3:ncol(compare88770)) {
  fit <- coxph(Surv(OS.time, OS) ~ compare88770[, i], data = compare88770)
  sum.surv <- summary(fit)
  c_index <- sum.surv$concordance
  ci <- confint(fit)
  lower_ci <- ci[1]
  upper_ci <- ci[2]
  p_value <- sum.surv$coefficients[,"Pr(>|z|)"]
  cindex88770 <- rbind(cindex88770, data.frame(variable = colnames(compare88770)[i], 
                                               cindex = c_index, 
                                               lower_ci = lower_ci, 
                                               upper_ci = upper_ci,
                                               p_value = p_value))
}
cindex88770 <- cindex88770 %>%
  group_by(variable) %>%
  summarise(
    C_index_sum = sum(cindex),
    CI_lower_avg = mean(lower_ci),
    CI_upper_avg = mean(upper_ci),
    p_value_avg = mean(p_value)
  )
save(cindex88770,file = 'cindex88770.RData')
##------筛选变量----------
comsam<-intersect(rownames(traincindex),rownames(cindex20685))
traincindex$variable<-rownames(traincindex)
traincindex<-as.data.frame(traincindex[comsam,])
cindex20685$variable<-rownames(cindex20685)
cindex20685<-as.data.frame(cindex20685[comsam,])
cindex42568$variable<-rownames(cindex42568)
cindex42568<-as.data.frame(cindex42568[comsam,])
cindex88770$variable<-rownames(cindex88770)
cindex88770<-as.data.frame(cindex88770[comsam,])

library(dplyr)
#--------------------绘制气泡图-----------------
##-----train-----
traincindex<-merge(traincindex,reference,by='variable')
traincindex<-traincindex[,-1]
colnames(traincindex)[6]<-'variable'
traincindex<- traincindex %>%
  mutate(significance = case_when(
    p_value_avg < 0.0001 ~ "****",
    p_value_avg < 0.001  ~ "***",
    p_value_avg < 0.01   ~ "**",
    p_value_avg < 0.05   ~ "*",
    TRUE            ~ ""  # pvalues 不显著时不显示标签
  ))
ggplot(traincindex, aes(x = `C_index_sum`, y = reorder(variable, `C_index_sum`), size = `C_index_sum`)) +
  geom_point(color = "#E64B35AA") +
  scale_size_continuous(range = c(2, 4)) +
  xlim(0.4, 1.0) +
  theme_classic() +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_text(aes(label = significance), color = "black", size = 2, hjust = -1) + # 设置 hjust 参数为 1.5
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
##---tcga-------
rownames(tcgacindex)<-tcgacindex$variable
tcgacindex<-tcgacindex[comsam,]
tcgacindex<-merge(tcgacindex,reference,by='variable')
tcgacindex<-tcgacindex[,-1]
colnames(tcgacindex)[6]<-'variable'
tcgacindex<- tcgacindex %>%
  mutate(significance = case_when(
    p_value_avg < 0.0001 ~ "****",
    p_value_avg < 0.001  ~ "***",
    p_value_avg < 0.01   ~ "**",
    p_value_avg < 0.05   ~ "*",
    TRUE            ~ ""  # pvalues 不显著时不显示标签
  ))
ggplot(tcgacindex, aes(x = `C_index_sum`, y = reorder(variable, `C_index_sum`), size = `C_index_sum`)) +
  geom_point(color = "#8491B4AA") +
  scale_size_continuous(range = c(2, 4)) +
  xlim(0.4, 1.0) +
  theme_classic() +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_text(aes(label = significance), color = "black", size = 2, hjust = -1) + # 设置 hjust 参数为 1.5
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
##---metacindex-----
rownames(metacindex)<-metacindex$variable
metacindex<-metacindex[comsam,]
metacindex<-merge(metacindex,reference,by='variable')
metacindex<-metacindex[,-1]
colnames(metacindex)[5]<-'variable'
metacindex<- metacindex %>%
  mutate(significance = case_when(
    p_value_avg < 0.0001 ~ "****",
    p_value_avg < 0.001  ~ "***",
    p_value_avg < 0.01   ~ "**",
    p_value_avg < 0.05   ~ "*",
    TRUE            ~ ""  # pvalues 不显著时不显示标签
  ))
ggplot(metacindex, aes(x = `C_index_sum`, y = reorder(variable, `C_index_sum`), size = `C_index_sum`)) +
  geom_point(color =  "#FFD54FAA") +
  scale_size_continuous(range = c(2, 4)) +
  xlim(0.4, 1.0) +
  theme_classic() +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_text(aes(label = significance), color = "black", size = 2, hjust = -1) + # 设置 hjust 参数为 1.5
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
##---20685-----
rownames(cindex20685)<-cindex20685$variable
cindex20685<-cindex20685[comsam,]
cindex20685<-merge(cindex20685,reference,by='variable')
cindex20685<-cindex20685[,-1]
colnames(cindex20685)[5]<-'variable'
cindex20685<- cindex20685 %>%
  mutate(significance = case_when(
    p_value_avg < 0.0001 ~ "****",
    p_value_avg < 0.001  ~ "***",
    p_value_avg < 0.01   ~ "**",
    p_value_avg < 0.05   ~ "*",
    TRUE            ~ ""  # pvalues 不显著时不显示标签
  ))
ggplot(cindex20685, aes(x = `C_index_sum`, y = reorder(variable, `C_index_sum`), size = `C_index_sum`)) +
  geom_point(color = '#00A087B2') +
  scale_size_continuous(range = c(2, 4)) +
  xlim(0.4, 1.0) +
  theme_classic() +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_text(aes(label = significance), color = "black", size = 2, hjust = -1) + # 设置 hjust 参数为 1.5
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
##---42568-----
rownames(cindex42568)<-cindex42568$variable
cindex42568<-cindex42568[comsam,]
cindex42568<-merge(cindex42568,reference,by='variable')
cindex42568<-cindex42568[,-1]
colnames(cindex42568)[6]<-'variable'
cindex42568<- cindex42568 %>%
  mutate(significance = case_when(
    p_value_avg < 0.0001 ~ "****",
    p_value_avg < 0.001  ~ "***",
    p_value_avg < 0.01   ~ "**",
    p_value_avg < 0.05   ~ "*",
    TRUE            ~ ""  # pvalues 不显著时不显示标签
  ))
ggplot(cindex42568, aes(x = `C_index_sum`, y = reorder(variable, `C_index_sum`), size = `C_index_sum`)) +
  geom_point(color = '#F39B7FB2') +
  scale_size_continuous(range = c(2, 4)) +
  xlim(0.4, 1.0) +
  theme_classic() +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_text(aes(label = significance), color = "black", size = 2, hjust = -1) + # 设置 hjust 参数为 1.5
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
##---88770----
rownames(cindex88770)<-cindex88770$variable
cindex88770<-cindex88770[comsam,]
cindex88770<-merge(cindex88770,reference,by='variable')
cindex88770<-cindex88770[,-1]
colnames(cindex88770)[5]<-'variable'
cindex88770<- cindex88770 %>%
  mutate(significance = case_when(
    p_value_avg < 0.0001 ~ "****",
    p_value_avg < 0.001  ~ "***",
    p_value_avg < 0.01   ~ "**",
    p_value_avg < 0.05   ~ "*",
    TRUE            ~ ""  # pvalues 不显著时不显示标签
  ))
ggplot(cindex88770, aes(x = `C_index_sum`, y = reorder(variable, `C_index_sum`), size = `C_index_sum`)) +
  geom_point(color = '#4DBBD5B2') +
  scale_size_continuous(range = c(2, 4)) +
  xlim(0.4, 1.0) +
  theme_classic() +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_text(aes(label = significance), color = "black", size = 2, hjust = -1) + # 设置 hjust 参数为 1.5
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
##----------48390-------
rownames(cindex48390)<-cindex48390$variable
cindex48390<-as.data.frame(cindex48390)
cindex48390<-cindex48390[comsam,]
cindex48390<-merge(cindex48390,reference,by='variable')
cindex48390<-cindex48390[,-1]
colnames(cindex48390)[5]<-'variable'
cindex48390<- cindex48390 %>%
  mutate(significance = case_when(
    p_value_avg < 0.0001 ~ "****",
    p_value_avg < 0.001  ~ "***",
    p_value_avg < 0.01   ~ "**",
    p_value_avg < 0.05   ~ "*",
    TRUE            ~ ""  # pvalues 不显著时不显示标签
  ))
ggplot(cindex48390, aes(x = `C_index_sum`, y = reorder(variable, `C_index_sum`), size = `C_index_sum`)) +
  geom_point(color = '#91D1C2FF') +
  scale_size_continuous(range = c(2, 4)) +
  xlim(0.4, 1.0) +
  theme_classic() +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_text(aes(label = significance), color = "black", size = 2, hjust = -1) + # 设置 hjust 参数为 1.5
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
