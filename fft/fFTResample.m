function x_t2 = fFTResample(x_t1,t1,t2)

resample x_t1 to x_t2 based on fFourierTransform, which manages NaNs
- create a full spectrum InversFT(t) of FT(x_t1)
- take x_t2 values of InversFT(t) at t2 sample points