# Calculation of site-specific evolutionary rates using the rosetta-ddg-based stability model
# and the WCN-based stress model. Comparison with empirical rates and calculation of 
# protein-by-protein goodness-of-fit (gof) measures
# OUTPUT:
# goodness-of-fit measures in "rosetta_gof.csv"
# empirical and predicted rates in "rosetta_rate_profiles.csv"

# load functions
source("misc_functions.r") # general helper functions 
source("model_functions.r") # rate predictors

# parameters
predictor <- "stability_model" # predictor function to use (from model_functions.r)
predict.rate = function(...) predict.rate.stability_model(...) # choose rate-predictor
print(c("predictor:",predictor))

# input file names
huang.fname <- "../Huang_et_al_data/wcn_profiles_209_monomers.csv" # wcn profiles
ej.fname = "../empirical_rate_calculations/evolutionary_rates/all_evolutionary_rates.csv" # empirical rates 

# output file names
profiles.fname <- "rosetta_rate_profiles.csv"
results.fname = "rosetta_gof.csv"

# input 
profiles.data = read.csv(huang.fname) # read sequence and structural profiles 
ej.data = read.csv(ej.fname) # read empirical rates
pdbs = substr(list.files("../ddG_calculations/rosetta/rosetta_ddG"),start=1,stop=4)



results = data.frame()
missing = data.frame()

for ( pdb.id in pdbs)
{
    print( pdb.id )

    chain = profiles.data$chain[profiles.data$pdb==pdb.id][1]
    ddG.data = read.ddG.rosetta( pdb.id, chain )


    if ( !is.null(ddG.data) )
    {
        site.index = (profiles.data$pdb==pdb.id)
        if ( nrow(ddG.data) != length(profiles.data$pdb[site.index]) )
        {
            print ("length mismatch")
            missing = rbind( missing, data.frame( pdb=pdb.id, issue="length mismatch") )
        }
        else
        {
            # use JC rates from ej.data
            r4s = ej.data[ej.data$pdb == pdb.id,"r4s_JC"]
            nr4s = r4s/mean(r4s)
            profiles.data$nr4s[site.index] = nr4s

            f = function(scale)
                    abs(cor(apply(ddG.data,
                        1,
                        function(x) predict.rate(x, scale=scale)),
                        nr4s,use="complete.obs"))
            start.time = Sys.time()
            o=optimize(f, interval=c(.2,10), maximum=T)
            end.time = Sys.time()
            elapsed.time = end.time - start.time
            print(c("elapsed.time",elapsed.time))
            scale.ddg=o$maximum
            cor.ddg=o$objective
            rate.ddg = apply(ddG.data, 1, function(x) predict.rate(x, scale=scale.ddg))
        
            method='p'
            profiles.data$ddg[site.index] = rate.ddg
            profiles.data$zddg[site.index] = scale(rate.ddg)

#   add rates predicted "locally" (i.e. for each protein) using ddg, lms, and ddg+lms
            lms  = profiles.data$wcn[site.index] #use the pfanm stress

            fit.ddgPlms <- lm(nr4s~rate.ddg+lms)
            fit.lms <- lm(nr4s~lms)
            fit.ddg <- lm(nr4s~rate.ddg)
     
# add fitted values to data.frame
            profiles.data[site.index,"rate.ddgPlms"] = fit.ddgPlms$coefficients[1]+
                                                  fit.ddgPlms$coefficients[2]*rate.ddg +
                                                  fit.ddgPlms$coefficients[3]*lms
            profiles.data[site.index,"rate.lms"] = fit.lms$coefficients[1]+fit.lms$coefficients[2]*lms
            profiles.data[site.index,"rate.ddg"] = fit.ddg$coefficients[1]+fit.ddg$coefficients[2]*rate.ddg

            profiles.data[site.index,"z.rate.lms"] = scale(profiles.data[site.index,"rate.lms"])
            profiles.data[site.index,"z.rate.ddg"] = scale(profiles.data[site.index,"rate.ddg"])
            profiles.data[site.index,"z.rate.ddgPlms"] = scale(profiles.data[site.index,"rate.ddgPlms"])

# goodness of fit measures 
# R2
            r2.ddgPlms <- summary(fit.ddgPlms)$r.squared
            r2.lms <- summary(fit.lms)$r.squared
            r2.ddg <- summary(fit.ddg)$r.squared
# Unique-Common analysis
            unique.ddg = r2.ddgPlms - r2.lms
            unique.lms = r2.ddgPlms - r2.ddg
            common.ddgORlms = r2.ddg - unique.ddg

# AICc per site (to compare between different proteins, see Roger 2014)
            nsites = length(nr4s)
            aic.ddgPlms <- AICc(fit.ddgPlms,ndat=nsites,kextra=1)/nsites # scale is 1 extra parameter
            aic.lms <- AICc(fit.lms,ndat=nsites,kextra=0)/nsites 
            aic.ddg <- AICc(fit.ddg,ndat=nsites,kextra=1)/nsites # scale  is 1 extra parameter

            row = data.frame( pdb=pdb.id,
                    scale.ddg=scale.ddg,
                    cor.ddg=cor.ddg,
                    cor.lms=cor(nr4s, profiles.data$rate.lms[site.index], method=method,use="complete.obs"),
                    cor.ddgPlms=cor(nr4s, profiles.data$rate.ddgPlms[site.index], method=method,use="complete.obs"),
                    r2.ddgPlms,r2.lms,r2.ddg,
                    unique.ddg,unique.lms,common.ddgORlms,
                    aic.ddgPlms,aic.lms,aic.ddg
                    )
            results = rbind(results, row)
            print( row )
        }
    }
    else
    {
        print ("ddG data missing")
        missing = rbind( missing, data.frame( pdb=pdb.id, issue="ddG data missing") )
    }
}

d=na.omit(profiles.data)

# output d and results to use for later analysis
write.csv(results,results.fname)
write.csv(d,profiles.fname)
