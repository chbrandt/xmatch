class Spline(object):
    '''
    Keep the state of a spline fit for a certain 'x' and 'y'
    '''
    _tck = None
    _q = 1

    def fit(self,x,y,normal=False):
        from scipy.interpolate import splrep
        tck = splrep(x, y, s=0)
        self._tck = tck
        if normal:
            self.normalize()

    @property
    def Q(self):
        return self._q

    @property
    def nodes(self):
        if self._tck is None:
            print("None curve has been fit. Try 'Spline.fit()' first.")
            return None
        return self._tck[0]

    @property
    def coefficients(self):
        if self._tck is None:
            print("None curve has been fit. Try 'Spline.fit()' first.")
            return None
        return self._tck[1]

    @property
    def degree(self):
        if self._tck is None:
            print("None curve has been fit. Try 'Spline.fit()' first.")
            return None
        return self._tck[2]

    def evaluate(self,x):
        if self._tck is None:
            print("None curve has been fit. Try 'Spline.fit()' first.")
            return None
        from scipy.interpolate import splev
        return splev(x,self._tck,der=0)

    def integrate(self,a,b):
        if self._tck is None:
            print("None curve has been fit. Try 'Spline.fit()' first.")
            return None
        from scipy.interpolate import splint
        return splint(a,b,self._tck)

    def normalize(self,a=None,b=None):
        if self._tck is None:
            print("None curve has been fit. Try 'Spline.fit()' first.")
            return None
        if a is None:
            a = self.nodes.min()
        if b is None:
            b = self.nodes.max()
        assert a<b, "Integration limits are inverted(?): a={};b={}".format(a,b)
        q = self.integrate(a,b)
        t,c,k = self._tck
        c = c/q
        self._tck = (t,c,k)
        self._q = q


def compile_interpolation_function(x,y,normal=False):
    '''
    Return Spline instance with 'x' and 'y' fit from a spline
    '''
    spl = Spline()
    spl.fit(x,y,normal=normal)
    return spl
# def compile_interpolation_function(x,y):
#     from scipy import interpolate
#     tck = interpolate.splrep(x, y, s=0)
#     def f(x,tck=tck):
#         return interpolate.splev(x, tck, der=0)
#     return f


def surface_density(data, surface_area, return_hist=False):
    from numpy import histogram,diff
    assert data.ndim == 1
    #print(data, surface_area)
    vec = data
    area_total = surface_area

    min_ = vec.min()
    max_ = vec.max()
    #print(min_,max_)
    #print(vec)

    h,b = histogram(vec, bins='auto', range=(min_,max_))
    h = h/area_total

    x = diff(b)/2 + b[:-1]
    func_nm = compile_interpolation_function(x,h)

    if return_hist:
        return func_nm,h,b
    return func_nm

def surface_densities(sample,columns,area_total,return_hist=True):
    '''
    Return pre-set density functions for each of 'columns'
    '''
    hist = {}
    bins = {}
    func_nm = {}
    for col in columns:
        vec = sample[col]
        func,h,b = surface_density(vec,area_total,True)
        func_nm[col] = func
        hist[col] = h
        bins[col] = b

    if return_hist:
        return func_nm,hist,bins
    return func_nm
