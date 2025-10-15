

# Read the data: Import on the right
attach(Zinc_m1)
# List of variables
names(Zinc_m1)
# Show first lines of data
head(Zinc_m1)
Zinc_m1[1:10,]
# Defining variables
price <- Price
d.price <- diff(price)
t <- Date
# Descriptive statistics
summary(price)
summary(d.price)
plot(price, type="o")
plot(d.price, type="o")

# Augmented Dickey-Fuller test (null hypothesis = data is non-stationary -->p<0.05: data is stationary)
adf.test(price, alternative="stationary")
adf.test(d.price, alternative="stationary")
# Kwiatkowski-Phillips-Schmidt-Shin (KPSS)  test (null hypothesis = data is stationary -->p<0.05: data is nonstationary)
kpss.test(price)

# ACF and PACF
acf(price)
pacf(price)

acf(d.price)
pacf(d.price)

# ARIMA(1,0,0) or AR(1)
arima(price, order=c(1,0,0))

# ARIMA(2,0,0) or AR(2)
arima(price, order=c(2,0,0))

# ARIMA(4,0,0) or AR(4)
arima(price, order=c(4,0,0))

# ARIMA(0,0,1) or MA(1)
arima(price, order=c(0,0,1))

# ARIMA(1,0,1) or AR(1)MA(1)
arima(price, order=c(1,0,1))

# ARIMA(1,1,0)
arima(d.price, order=c(1,0,0))

# ARIMA(0,1,1)
arima(d.price, order=c(0,0,1))

# ARIMA(1,1,1)
arima(d.price, order=c(1,0,1))

# ARIMA(1,1,3)
arima(d.price,order=c(1,0,3))

# ARIMA(2,1,3)
arima(d.price,order=c(2,0,3))




