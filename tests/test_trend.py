import unittest
import trend
import numpy.ma as ma
import numpy as np


class TestWindow (unittest.TestCase) : 
    def setUp(self) : 
        # global, 0.5 deg dataset
        self.lats = np.arange(-89.75, 90, 0.5) 
        self.lons = np.arange(-179.75, 180, 0.5)
        
        self.geog_box = [ -9.25, 179.75 , 20.25, 74.75 ] 
        self.lat_size_win = int(((74.75-20.25)*2)+1)
        self.lon_size_win = int(((179.75-(-9.25))*2)+1)
        self.w = trend.Window(self.lats, self.lons, self.geog_box) 
        
    def test_dataset_size(self) : 
        """global 0.5deg dataset should be 720x360"""
        self.assertEqual(self.w._dataset_shape[0], 360.)
        self.assertEqual(self.w._dataset_shape[1], 720.)
        
    def test_window_size(self) : 
        self.assertEqual(self.w._window_shape[0], self.lat_size_win)
        self.assertEqual(self.w._window_shape[1], self.lon_size_win)
        
    def test_window_geom(self) : 
        """ensures window located correctly in dataset"""
        lat_slice = self.w._window[0]
        lon_slice = self.w._window[1]
        self.assertEqual(self.lats[lat_slice.start], 20.25)
        self.assertEqual(self.lats[lat_slice.stop-1], 74.75)
        self.assertEqual(self.lons[lon_slice.start], -9.25)
        self.assertEqual(self.lons[lon_slice.stop-1], 179.75)
        
    def test_set_window(self) : 
        window_data = np.ones( (self.lat_size_win, self.lon_size_win) )
        x = self.w.set_window(window_data)
        
        # check output geometry
        self.assertEqual(x.shape[0], 360)
        self.assertEqual(x.shape[1], 720)
        
        # check output is masked
        self.assertTrue(ma.is_masked(x))
        
        # check that the window is only thing in returned array
        win_masked = ma.count_masked(x)
        win = ma.count(x)
        self.assertEqual(win, window_data.size)
        self.assertEqual(win_masked, x.size - window_data.size)
        self.assertTrue(np.all(x[self.w._window] == window_data))
        
    def test_get_window(self) : 
        dataset = np.zeros( (360,720) ) 
        dataset[self.w._window] = 1
        
        x = self.w.get_window(dataset)
        
        # check window geometry
        self.assertEqual(x.shape[0], self.lat_size_win)
        self.assertEqual(x.shape[1], self.lon_size_win)
        
        # check no zeros in returned data
        self.assertEqual(np.count_nonzero(x==0),0)
        
        
        
    
        
        