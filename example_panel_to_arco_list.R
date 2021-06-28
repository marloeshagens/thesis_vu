# = Generate a small panel as example = #
set.seed(123)
time=sort(rep(1:100,2))
unit=rep(c("u1","u2"),100)
v1=rnorm(200)
v2=rnorm(200)
panel=data.frame(time=time,unit=unit,v1=v1,v2=v2)
head(panel)

data=panel_to_ArCo_list(panel,time="TIME_PERIOD",unit="geo",variables = "siec")
head(data$v1)
