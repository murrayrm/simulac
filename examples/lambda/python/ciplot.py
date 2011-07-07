def ciplot(x, lower, upper, **keywords):

    import numpy as np
    import matplotlib.pyplot as mp
     
    # ciplot(lower,upper)       
    # ciplot(lower,upper,x)
    # ciplot(lower,upper,x,colour)
    #
    # Plots a shaded region on a graph between specified lower and upper
    # confidence intervals (L and U).  l and u must be vectors of the same
    # length.  Uses the 'fill' function, not 'area'. Therefore multiple
    # shaded plots can be overlayed without a problem. Make them
    # transparent for total visibility.
    #
    # Use keywords to set color of plot.
    #
    # Modified from ciplot.m, by Raymond Reynolds 24/11/06

    if (lower.shape[0] != upper.shape[0]):
        fprintf(stderr, 'lower and upper vectors must be same length');
    else:
        mp.fill(
            np.append(x, x[::-1]),
            np.append(upper, lower[::-1]), **keywords)
