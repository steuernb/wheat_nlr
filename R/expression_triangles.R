
library(ggtern)

 f <- read.delim("steuernb/nlr_paper/v1.1_analysis/rnaseq/dge/compare_F.txt",row.names="gene")
 c <- read.delim("steuernb/nlr_paper/v1.1_analysis/rnaseq/dge/compare_C.txt",row.names="gene")
 h <- read.delim("steuernb/nlr_paper/v1.1_analysis/rnaseq/dge/compare_H.txt",row.names="gene")
 hcf <- read.delim("steuernb/nlr_paper/v1.1_analysis/rnaseq/dge/compare_HCF.txt",row.names="gene")


#############################
##  Flagellin

x <- f$H000[order(f$F030)]
y <- f$F030[order(f$F030)]
z <- f$F180[order(f$F030)]
t<- f$maxTPM[order(f$F030)]

df<- data.frame(x=x, y=y, z=z, t=t)

p1<-ggtern(data=df,aes(x,y,z)) +
  geom_point(aes(color = t,size = t  ) ) +
  scale_size_continuous(range = c(0.5, max(df$t)/10)) +
  scale_color_gradient(low = 'yellow', high = 'red', limits=c(0,67)) +
  theme_bw() + 
  theme(legend.position="none") +
  theme_nogrid() +
  labs(x="H000",y="F030",z="F180",title="Flagellin") 


#############################
##  Chitin

x <- c$H000[order(c$C030)]
y <- c$C030[order(c$C030)]
z <- c$C180[order(c$C030)]
t<- c$maxTPM[order(c$C030)]
mysize <- round(t) 
mycol <- round(t/5)*5
df<- data.frame(x=x, y=y, z=z, t=t)
p2<-ggtern(data=df,aes(x,y,z)) +
  geom_point(aes(color = t,size = t  ) ) +
  scale_size_continuous(range = c(0.5, max(df$t)/10)) +
  scale_color_gradient(low = 'yellow', high = 'red', limits=c(0,67)) +
  theme_bw() + 
  theme(legend.position="none") +
  theme_nogrid() +
  labs(x="H000",y="C030",z="C180",title="Chitin") 


############################
## Mock


x <- h$H000[order(h$H030)]
y <- h$H030[order(h$H030)]
z <- h$H180[order(h$H030)]
t<-  h$maxTPM[order(h$H030)]
mysize <- round(t) 
mycol <- round(t/5)*5
df<- data.frame(x=x, y=y, z=z, t=t)
p3<-ggtern(data=df,aes(x,y,z)) +
  geom_point(aes(color = t,size = t  ) ) +
  scale_size_continuous(range = c(0.5, max(df$t)/10)) +
  scale_color_gradient(low = 'yellow', high = 'red', limits=c(0,67)) +
  theme_bw() +
  theme(legend.position="none") +
  theme_nogrid() +
  labs(x="H000",y="H030",z="H180",title="Mock") 

##########################
x <- hcf$H030[order(hcf$F030)]
y <- hcf$C030[order(hcf$F030)]
z <- hcf$F030[order(hcf$F030)]
t<- hcf$maxTPM[order(hcf$F030)]
mysize <- round(t) 
mycol <- round(t/5)*5
df<- data.frame(x=x, y=y, z=z, t=t)
p4<-ggtern(data=df,aes(x,y,z)) +
  geom_point(aes(color = t,size = t  ) ) +
  scale_size_continuous(range = c(0.5, max(df$t)/10)) +
  scale_color_gradient(low = 'yellow', high = 'red', limits=c(0,67)) +
  theme_bw() +
  theme(legend.position="none") +
  theme_nogrid() +
  labs(x="H030",y="C030",z="F030",title="Treatment") 




#svg(file="steuernb/nlr_paper/v1.1_analysis/rnaseq/dge/triangles.svg",  width=6.7, height=3.5)
png(file="steuernb/nlr_paper/v1.1_analysis/rnaseq/dge/triangles.png", width=40, height=20,units="cm", res=300 )

grid.arrange(p1, p2,p3,p4, nrow = 2)
dev.off()
