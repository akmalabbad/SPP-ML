============== COVARIATES ===============
# From RDA Files
Depth : The depth of the center of earthquake
Dip   : Dip (The angle formed by the fault plane with the horizontal plane and
measured in the vertical plane)
Ses   : Earthquake distance to nearest fault (100 km)
Sub   : Earthquake distance to nearest subduction zones (100 km)
Vul   : Earthquake distance to nearest volcanoes (100 km)
Strike: Strike (The direction of the horizontal line that lies in the fault plane is
measured clockwise from the north)

# From CSV File
Long  : Longitude
Lat   : Latitude

# Code Explanation
Vol               : The area size divided total number of dummy points. If it is rectangle. Otherwise, it can be generated automatically using `quadscheme`
F_prime           : This is the Ripley’s K-function in equation 5 from `XGBoostPP: Tree-based Estimation of Point Process
Intensity Functions` paper
dynamic_weighted  : Our poisson loss can accept poisson loss and weighted poisson loss 


# TODO:
- Buat formula vol dan F_prime di supplementary.R
- Default vol = 0 pada LGBPP agar weight = 1
- Buat if else di Gempa Tuning untuk accomodate dynamic_weighted = True
- Buat parameter loss option either poisson or logistic, jadi langkah di atas tidak perlu yang dynamic_weighted
- Tambahkan logistic ke LGBPP
- Jangan lupa quadscheme dan quadscheme.logi di if else
- Bandingkan keempat performanya
- Map nya belum pas di pulau sumatera



