# script that is going to make all the final figures in the paper

# input and output folders
input.dir <- "."
output.dir <- "."

# required libraries
require("ggplot2")
require("MASS")


# foldx vs. rosetta
# still need to investigate warnings in "one-param-model-rosetta.r", can
# they be ignored?
if (T)
{
    # results by protein
    MC_foldx = read.csv(file.path(input.dir,"foldx_gof.csv"))
    MC_rosetta = read.csv(file.path(input.dir,"rosetta_gof.csv"))

    pdf(file.path(output.dir,"foldx-vs-rosetta.pdf"), width=7, height=6)
    par(mar=c(5, 5, 3, 1))

    i = match(MC_rosetta$pdb, MC_foldx$pdb)
    plot(MC_rosetta$cor.ddg, MC_foldx[i,]$cor.ddg,
        xlab='correlation, Rosetta',
        ylab='correlation, FoldX',
        xlim=c(.5, .75),
        ylim=c(.5, .75),
        pch=19,
        col='blue',
        axes=F,
        cex.lab=1.2)
    abline(0, 1)
    axis(1, cex.axis=1.2)
    axis(2, cex.axis=1.2)
    text(MC_rosetta$cor.ddg-.008, MC_foldx[i,]$cor.ddg+.007, MC_rosetta$pdb)
    dev.off()
    
    cat("Shift between foldx model and rosetta:\n")
    print( t.test( MC_foldx[i,]$cor.ddg, MC_rosetta$cor.ddg, paired=T ) )
}

# correlations by protein
if (T)
{
    # results by protein
    MC_foldx = read.csv(file.path(input.dir,"foldx_gof.csv"))
    MC_foldx_pf = read.csv(file.path(input.dir,"pf_foldx_gof.csv"))
    
    pdf(file.path(output.dir,"rate-correlations.pdf"), width=14, height=6)
    split.screen(c(1, 2))

    screen(1)
    par(mar=c(5, 5, 3, 1))

    plot(MC_foldx$scale.ddg, MC_foldx$cor.ddg,
        xlab=expression(paste('scale ', alpha)),
        ylab=expression(paste('correlation ', italic(r), ', ', alpha, ' fitted')),
        xlim=c(0, 3.5),
        ylim=c(.2, .75),
        pch=19,
        col='blue',
        axes=F,
        cex.lab=1.5
    )
    axis(1, cex.axis=1.5)
    axis(2, cex.axis=1.5)
    mtext("A", side=3, at=-.5, line=-1, cex=1.5, font=2)
    
    cat("Correlation between scale.ddg and cor.ddg:\n")
    print( cor.test( MC_foldx$scale.ddg, MC_foldx$cor.ddg ) )
    
    screen(2)
    par(mar=c(5, 5, 3, 1))

    plot(MC_foldx_pf$cor.ddg, MC_foldx$cor.ddg,
        xlab=expression(paste('correlation ', italic(r), ', ', alpha, ' = 1')),
        ylab=expression(paste('correlation ', italic(r), ', ', alpha, ' fitted')),
        xlim=c(.2, .75),
        ylim=c(.2, .75),
        pch=19,
        col='blue',
        axes=F,
        cex.lab=1.5
    )
    axis(1, cex.axis=1.5)
    axis(2, cex.axis=1.5)
    abline(0, 1)
    mtext("B", side=3, at=.1, line=-1, cex=1.5, font=2)
    
    cat("Shift between MC model and MC_pf model:\n")
    print( t.test( MC_foldx$cor.ddg, MC_foldx_pf$cor.ddg, paired=T ) )
    
    close.screen(all.screens=T)
    dev.off()
}


# ddg vs. lms
if (T)
{
    # results by protein
    MC_foldx = read.csv(file.path(input.dir,"foldx_gof.csv"))
    
    pdf(file.path(output.dir,"ddg-vs-lms.pdf"), width=7, height=6)
    
    par(mar=c(5, 5, 3, 1))

    plot(MC_foldx$cor.lms, MC_foldx$cor.ddg,
        #xlab='correlation, lms',
        xlab=expression(paste('correlation, ', Delta, Delta, G^"*")),
        ylab=expression(paste('correlation, ', Delta, Delta, G)),
        xlim=c(.2, .75),
        ylim=c(.2, .75),
        pch=19,
        col='blue',
        axes=F,
        cex.lab=1.2
    )
    axis(1, cex.axis=1.2)
    axis(2, cex.axis=1.2)
    abline(0, 1)
    
    dev.off()
    
    cat("Correlation between lms model and ddg model:\n")
    print( cor.test( -MC_foldx$cor.lms, MC_foldx$cor.ddg, paired=T ) )
    
    cat("Shift between lms model and ddg model:\n")
    print( t.test( -MC_foldx$cor.lms, MC_foldx$cor.ddg, paired=T ) )
}


# rate distribution plot
if (F)
{
    d.fname = file.path(input.dir, paste("foldx_rate_profiles.csv", sep=""))
    d = read.csv(d.fname)
    
    # calculate loess regression
    r4s.loess.fun = approxfun(lowess(d$rate.ddg, d$nr4s))
    x = seq( -1, 3, .05)
    r4s.loess = data.frame(x=x, y=r4s.loess.fun(x))
    
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
            geom_abline(intercept = 0, aes(color='x=y'), slope = 1) + 
            # loess approx. of r4s vs rate.ddg
            geom_line(data=r4s.loess, aes(x=x, y=y, color='LOESS fit'), size=1.5) +
            # legend
            scale_color_manual(values=c("x=y"="black",
                                        "LOESS fit"="red")) +
            guides(color=guide_legend(title=NULL,
                        override.aes=list(size=c(1.5, .5)))
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
                  legend.text=element_text(size=14)
                  ) +
            # axis ticks and labels
            scale_x_continuous(breaks=c(0, .5, 1, 1.5, 2), 
                               labels=c("0.0", "0.5", "1.0", "1.5", "2.0")) +
            scale_y_continuous(breaks=c(0, .5, 1, 1.5, 2, 2.5), 
                               labels=c("0.0", "0.5", "1.0", "1.5", "2.0", "2.5")) +
            xlab(expression(paste("predicted rate (", Delta, Delta, italic(G), ")"))) +
            ylab("measured rate (Rate4Site)")
    fname = file.path(output.dir, "MC_rates_distribution.pdf")
    ggsave(plot=plot,
        fname,
        width=7.5,
        height=5)
}
