"""Create grid inside Polygon objects."""

import geopandas as gpd
import pandas as pd
import numpy as np
import math
from shapely.geometry import Point

from sklearn.neighbors.kde import KernelDensity
from scipy.spatial.distance import cdist


class Grid(object):
    """
    Construction process of a Grid object
    """

    def __init__(self, data_type, step, bandwidth):

        # Polygon or Point
        self.data_type = data_type 

        # the smaller the step, the denser the grid, the higher the accuracy
        self.step = step
        
        # Bandwidth used for kernel density estimation
        self.bandwidth = bandwidth

    # Editted from https://github.com/lgervasoni/urbansprawl
    def create_grid_from_gdf(self, object_gdf):
        """
        Create a regular grid with the input geodataframe
        Editted from https://github.com/lgervasoni/urbansprawl

        Returns
        ----------
        geopandas.GeoDataFrame
        """
        # Get bounding box
        west, south, east, north = object_gdf.total_bounds
        # Create indices
        grid_gdf = gpd.GeoDataFrame( [ Point(i,j) for i in np.arange(west, east, self.step) for j in np.arange(south, north, self.step) ], columns=["geometry"] )
        # Set projection
        grid_gdf.crs = object_gdf.crs

        return grid_gdf    


    def WeightedKernelDensityEstimation(self, X, Weights, Y, max_mb_per_chunk = 1000):
        """ 
        Computes a Weighted Kernel Density Estimation
        Editted from https://github.com/lgervasoni/urbansprawl
        
        Returns
        ----------
        pd.Series
            returns an array of the estimated densities rescaled between [0;1]
        """
        def get_megabytes_pairwise_distances_allocation(X, Y):
            # Calculate MB needed to allocate pairwise distances
            return len(X) * len(Y) * 1e-6
        
        # During this procedure, pairwise euclidean distances are computed between inputs points X and points to estimate Y
        # For this reason, Y is divided in chunks to avoid big memory allocations. At most, X megabytes per chunk are allocated for pairwise distances
        Y_split = np.array_split( Y, math.ceil( get_megabytes_pairwise_distances_allocation(X,Y) / max_mb_per_chunk ) )
        
        """
        ### Step by step
        # Weighed KDE: Sum{ Weight_i * K( (X-Xi) / h) }
        W_norm = np.array( Weights / np.sum(Weights) )
        cdist_values = cdist( Y, X, 'euclidean') / bandwidth
        Ks = np.exp( -.5 * ( cdist_values ) ** 2  )
        PDF = np.sum( Ks * W_norm, axis=1)
        """
        """
        ### Complete version. Memory consuming
        PDF = np.sum( np.exp( -.5 * ( cdist( Y, X, 'euclidean') / bandwidth ) ** 2  ) * ( np.array( Weights / np.sum(Weights) ) ), axis=1)
        """

        ### Divide Y in chunks to avoid big memory allocations
        PDF = np.concatenate( [ np.sum( np.exp( -.5 * ( cdist( Y_i, X, 'euclidean') / self.bandwidth ) ** 2  ) * ( np.array( Weights / np.sum(Weights) ) ), axis=1) for Y_i in Y_split ] )
        
        # np.concatenate 数组拼接 
        # np.concatenate((a, b), axis=None)  array([1, 2, 3, 4, 5, 6])
            
        # Rescale
        return pd.Series( PDF / PDF.sum() )


    def calculate_kde(self, object_gdf, X_weights= None):
        """
        Computes a Kernel Density Estimation
        Editted from https://github.com/lgervasoni/urbansprawl

        Returns
        ----------
        pandas.Series
        """

        self.grid_gdf = self.create_grid_from_gdf( object_gdf = object_gdf)

        if self.data_type == "Polygon":
            X_b = [ [p.x,p.y] for p in object_gdf.geometry.centroid.values]

        if self.data_type == "Point":
            X_b = [ [p.x,p.y] for p in object_gdf.geometry.values]
        
        X = np.array(X_b)
        # Points where the probability density function will be evaluated
        Y = np.array( [ [p.x,p.y] for p in self.grid_gdf.geometry.values ] )
        
        if X_weights is None:

            kde = KernelDensity(kernel='gaussian', bandwidth = self.bandwidth).fit(X)

            # Sklearn returns the results in the form log(density)
            p = np.exp(kde.score_samples(Y))

        else:
        # Weighted Kernel Density Estimation   
            p = self.WeightedKernelDensityEstimation(X = X, Weights = X_weights, Y = Y,max_mb_per_chunk = 1000)

            return pd.Series( p / p.max() )

        # Returns the spatial distribution density
        return pd.Series( p / p.max() )

    
    def computed_grid(self, object_gdf, X_weights = None):
        """ 
        Calculate land use mix indices on input grid
        Editted from https://github.com/lgervasoni/urbansprawl

        Returns
        ----------
        pandas.DataFrame
        """
        self.grid_gdf = self.create_grid_from_gdf( object_gdf)

        if X_weights is None:
            self.grid_gdf["no_weighted"] = self.calculate_kde(object_gdf = object_gdf,
                                                            X_weights = None)
        else:
            self.grid_gdf["weighted"] = self.calculate_kde(object_gdf = object_gdf,
                                                            X_weights = X_weights)

        return self.grid_gdf


