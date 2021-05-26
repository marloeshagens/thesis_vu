

# TRY OUT - what happens with more variables?
data = data.q2
treated.unit = 1
t0 = 51
lag=0

for (i in 1:length(data)) {
  aux = length(unique(colnames(data[[i]])))
  k = ncol(data[[i]])
  if (aux < k) {
    colnames(data[[i]]) = paste("unit", 1:ncol(data[[i]]), 
                                sep = "")
  }
}


Y = Reduce("cbind", lapply(data, function(x) x[, treated.unit])) 
# takes first column (treated variable) of both variables and puts them next to each other; 
# now, Y is a matrix

X = Reduce("cbind", lapply(data, function(x) x[, -treated.unit]))
# takes all the rest of the columns (untreated units) of both variables and puts them next to each other;
# thus, we have var1: unit 2,3,4,5,6;var2: unit 2,3,4,5,6 in columns next to each other

aux = list()
for (i in 1:length(data)) {
  aux[[i]] = paste(names(data)[i], colnames(data[[i]])[-treated.unit], 
                   sep = ".")
}

colnames(X) = unlist(aux)
X

Y.raw = Y


T=nrow(X)
y.fit = matrix(Y[1:(t0 - 1 - lag),], ncol = length(data)) 
y.pred = matrix(Y[-c(1:(t0 - 1 - lag)), ], ncol = length(data))
x.fit = X[1:(t0 - 1 - lag), ]
x.pred = X[-c(1:(t0 - 1 - lag)), ]
save.cf = matrix(NA, nrow(y.pred), length(data))
save.fitted = matrix(NA, nrow(Y), length(data))
model.list = list()
for (i in 1:length(data)) {
  model = fn(x.fit, y.fit[, i],...)
  model.list[[i]] = model
  contra.fact = p.fn(model, x.pred)
  save.cf[, i] = contra.fact
  save.fitted[, i] = p.fn(model, X)
}