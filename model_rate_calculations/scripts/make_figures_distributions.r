# script that makes the measured-predicted comparisons

# input and output folders
input.dir <- "."
output.dir = "."

# input and output files
d.fname = file.path(input.dir, paste("foldx_rate_profiles.csv"))
ddg.fname = file.path(output.dir, "foldx_ddg_rates_distribution.pdf")
lms.fname = file.path(output.dir, "lms_rates_distribution.pdf")
ddgPlms.fname = file.path(output.dir, "foldx_ddgPlms_rates_distribution.pdf")

# parameters
read.d = T
plot.ddg=T
plot.lms=F
plot.ddgPlms=F

nbins = 20 # number of x-variable bins for nr4s averages
bin.method = "quantile" # "quantile" or "equal"
z99 = 2.576
zr = z99

# set colors
color.average = "#FFFF99"
color.fit = "black"

# required libraries
require("ggplot2")
require("MASS")
source("bin_functions.r")

# read d file
if (read.d) {
    d = read.csv(d.fname)
}


# rate distribution plot ddg
if (plot.ddg)
{
    # bin x-coordinate
    df = d[,c("rate.ddg","nr4s")]
    df2 = bin_df(df,col=1,nbins=nbins,method=bin.method)
    
    # calculate 2d density
    dens = MASS::kde2d(d$rate.ddg, d$nr4s, n = 100)
    densf = data.frame(expand.grid(x=dens$x, y=dens$y), z=as.vector(dens$z))
    
    plot = ggplot() +
        # plot 2d density of r4s and rate.ddg data
        geom_polygon(data=densf,
                     aes(x=x, y=y, z=z, fill=..level..),
                     stat="contour",
                     breaks=seq(min(densf$z), max(densf$z), length.out=100)) +
        scale_fill_gradient(low="aquamarine", high="midnightblue") +
        labs(fill="density") +
        
        # ddg-rate prediction is simple x=y line:
        geom_abline(aes(color="x=y"), intercept = 0, slope = 1) + 
        
        # bin-means +- sd
        geom_pointrange(data=df2, 
                        aes(x=mean.rate.ddg, 
                            y=mean.nr4s, 
                            ymin=p25.nr4s, 
                            ymax=p75.nr4s,
                            color='average'), 
                        size=.7) +
        # legend
        scale_colour_manual(values=c("x=y"=color.fit, "average"=color.average)) +
        #scale_shape_manual(values=c("x=y"=19, "average"=NA)) +
        guides(color=guide_legend(order=1,title=NULL,
                                  override.aes=list(size=c(.7, .7),shape=c(19,NA),background=element_rect())),
               fill=guide_legend(order=2)
        ) + 
        # x and y limits
        coord_cartesian(xlim = c(-.1, 2.1), ylim = c(-.1, 2.8) ) +
        # plot theme
        theme_classic() + # classic theme
        theme(axis.text.x=element_text(size=16),
              axis.text.y=element_text(size=16),
              axis.title.x=element_text(size=16, vjust=-.5),
              axis.title.y=element_text(size=16, vjust=.3),
              legend.title=element_text(size=14, face='plain'),
              legend.key=element_blank(),
              legend.text=element_text(size=14),
              legend.background=element_rect(fill="gray")
        ) +
        # axis ticks and labels
        scale_x_continuous(breaks=c(0, .5, 1, 1.5, 2), 
                           labels=c("0.0", "0.5", "1.0", "1.5", "2.0")) +
        scale_y_continuous(breaks=c(0, .5, 1, 1.5, 2, 2.5), 
                           labels=c("0.0", "0.5", "1.0", "1.5", "2.0", "2.5")) +
        xlab(expression(paste("predicted rate (", Delta, Delta, italic(G), ")"))) +
        ylab("measured rate (Rate4Site)")
    #debug
    ggsave(plot=plot, file=ddg.fname, width=7.5, height=5)
    print(plot)
}


# rate distribution plot lms
if (plot.lms)
{
    # bin x-coordinate
    df = d[,c("rate.lms","nr4s")]
    df2 = bin_df(df,col=1,nbins=nbins,method=bin.method)
    
    # calculate 2d density
    dens = MASS::kde2d(d$rate.lms, d$nr4s, n = 100)
    densf = data.frame(expand.grid(x=dens$x, y=dens$y), z=as.vector(dens$z))

    plot = ggplot() +
            # plot 2d density of r4s and rate.lms data
            geom_polygon(data=densf,
                aes(x=x, y=y, z=z, fill=..level..),
                stat="contour",
                breaks=seq(min(densf$z), max(densf$z), length.out=100)) +
            scale_fill_gradient(low="aquamarine", high="midnightblue") +
            labs(fill="density") +
        
            # lms-rate prediction is simple x=y line:
            geom_abline(aes(color="x=y"), intercept = 0, slope = 1) + 
        
            # bin-means +- sd
            geom_pointrange(data=df2, 
                            aes(x=mean.rate.lms, 
                                y=mean.nr4s, 
                                ymin=p25.nr4s, 
                                ymax=p75.nr4s,
                                color='average'), 
                            size=.7) +
            # legend
            scale_colour_manual(values=c("x=y"=color.fit, "average"=color.average)) +
            #scale_shape_manual(values=c("x=y"=19, "average"=NA)) +
            guides(color=guide_legend(order=1,title=NULL,
                        override.aes=list(size=c(.7, .7),shape=c(19,NA))),
                   fill=guide_legend(order=2)
                   ) + 
            # x and y limits
            coord_cartesian(xlim = c(-.1, 2.1), ylim = c(-.1, 2.8) ) +
            # plot theme
            theme_classic() + # classic theme
            theme(axis.text.x=element_text(size=16),
                  axis.text.y=element_text(size=16),
                  axis.title.x=element_text(size=16, vjust=-.5),
                  axis.title.y=element_text(size=16, vjust=.3),
                  legend.title=element_text(size=14, face='plain'),
                  legend.key=element_blank(),
                  legend.text=element_text(size=14),
                  legend.background=element_rect(fill="gray")
                  ) +
            # axis ticks and labels
            scale_x_continuous(breaks=c(0, .5, 1, 1.5, 2), 
                               labels=c("0.0", "0.5", "1.0", "1.5", "2.0")) +
            scale_y_continuous(breaks=c(0, .5, 1, 1.5, 2, 2.5), 
                               labels=c("0.0", "0.5", "1.0", "1.5", "2.0", "2.5")) +
            xlab(expression(paste("predicted rate (", italic("MLmS"),")"))) +
            ylab("measured rate (Rate4Site)")
#debug
    ggsave(plot=plot, file=lms.fname, width=7.5, height=5)
    print(plot)
}


# rate distribution plot ddgPlms
if (plot.ddgPlms)
{
    # bin x-coordinate
    df = d[,c("rate.ddgPlms","nr4s")]
    df2 = bin_df(df,col=1,nbins=nbins,method=bin.method)
    
    # calculate 2d density
    dens = MASS::kde2d(d$rate.ddgPlms, d$nr4s, n = 100)
    densf = data.frame(expand.grid(x=dens$x, y=dens$y), z=as.vector(dens$z))

    plot = ggplot() +
            # plot 2d density of r4s and rate.ddgPlms data
            geom_polygon(data=densf,
                aes(x=x, y=y, z=z, fill=..level..),
                stat="contour",
                breaks=seq(min(densf$z), max(densf$z), length.out=100)) +
            scale_fill_gradient(low="aquamarine", high="midnightblue") +
            labs(fill="density") +
        
            # ddgPlms-rate prediction is simple x=y line:
            geom_abline(aes(color="x=y"), intercept = 0, slope = 1) + 
        
            # bin-means +- sd
            geom_pointrange(data=df2, 
                            aes(x=mean.rate.ddgPlms, 
                                y=mean.nr4s, 
                                ymin=p25.nr4s, 
                                ymax=p75.nr4s,
                                color='average'), 
                            size=.7) +
            # legend
            scale_colour_manual(values=c("x=y"=color.fit, "average"=color.average)) +
            #scale_shape_manual(values=c("x=y"=19, "average"=NA)) +
            guides(color=guide_legend(title=NULL,override.aes=list(size=c(.7, .7),shape=c(19,NA)),order=1),
                   fill =guide_legend(order=2)
                   ) + 
            # x and y limits
            coord_cartesian(xlim = c(-.1, 2.1), ylim = c(-.1, 2.8) ) +
            # plot theme
            theme_classic() + # classic theme
            theme(axis.text.x=element_text(size=16),
                  axis.text.y=element_text(size=16),
                  axis.title.x=element_text(size=16, vjust=-.5),
                  axis.title.y=element_text(size=16, vjust=.3),
                  legend.title=element_text(size=14, face='plain'),
                  legend.key=element_blank(),
                  legend.text=element_text(size=14),
                  legend.background=element_rect(fill="gray")
                  ) +
            # axis ticks and labels
            scale_x_continuous(breaks=c(0, .5, 1, 1.5, 2), 
                               labels=c("0.0", "0.5", "1.0", "1.5", "2.0")) +
            scale_y_continuous(breaks=c(0, .5, 1, 1.5, 2, 2.5), 
                               labels=c("0.0", "0.5", "1.0", "1.5", "2.0", "2.5")) +
            xlab(expression(paste("predicted rate (", Delta, Delta, italic(G),"+",italic("MLmS"), ")"))) +
            ylab("measured rate (Rate4Site)")
#debug
    ggsave(plot=plot, file=ddgPlms.fname, width=7.5, height=5)
    print(plot)
}


