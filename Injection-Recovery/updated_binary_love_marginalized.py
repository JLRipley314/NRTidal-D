import numpy as np

def convert_lambda_s_to_lambda_a_marginalized(lambda_s,q):
    """
    Marginalized Binary Love relations. 
    
    Lambda_s = (Lambda_1 + Lambda_2)/2
    Lambda_a = (Lambda_1 - Lambda_2)/2
   
    The Binary Love relations give us Lambda_a(Lambda_s).

    This relation has some error \delta. 
    Following [1], we marginalize over this error by sampling from 
    a Gaussian distributions with mean/variance
    equal to the mean/variance of the distribution of predictions from the 
    binary Love relation for different EOS.
    Thus this function will give you a different answer if you call it multiple
    times with the same input, unless you manually reset the random seed every
    time you call this function.

    [1] https://arxiv.org/abs/1903.03909

    DOI: https://doi.org/10.1103/PhysRevD.99.083016
    """
    #define bmatrix and cmatrix
    bmatrix=np.array([[-14.4,14.45],[31.36,-32.25],[-22.44,20.35]])
    cmatrix=np.array([[-15.25,15.37],[37.33,-43.20],[-29.93,35.18]])
    
    #define other quantities for computing the binary love relation
    n=0.743
    alpha=1
    Fn=(1-q**(10/(3-n)))/(1+q**(10/(3-n)))
    #compute the sum in the numerator and denominator of binary love relation
    num_sum=0
    den_sum=0
    for i in range(3):
        for j in range(2):
            num_sum+=bmatrix[i,j]*q**(j+1)*lambda_s**(-(i+1)/5.)
            den_sum+=cmatrix[i,j]*q**(j+1)*lambda_s**(-(i+1)/5.)
    #compute binary love relation
    lambda_a=Fn*(1+num_sum)/(1+den_sum)*lambda_s**(alpha)

    #define coefficients mu's and sigma's
    mu_1=3.509e-3
    mu_2=9.351e-1
    mu_3=-18.07
    mu_4=27.56
    mu_5=-10.10
    sigma_1=-2.074e-7
    sigma_2=1.492e-3
    sigma_3=-4.891e-2
    sigma_4=8.207e-1
    sigma_5=-1.308
    sigma_6=-63.76
    sigma_7=11.14
    sigma_8=75.25
    sigma_9=-23.69 

    #construct the marginalized mean and standard deviation 
    mu_r_lambda_s=mu_1*lambda_s+mu_2
    mu_r_q=mu_3*q**2+mu_4*q+mu_5
    sigma_r_lambda_s=(sigma_1*lambda_s**(5./2)+sigma_2*lambda_s**(3./2) +sigma_3*lambda_s+sigma_4*lambda_s**(1./2)+sigma_5)
    sigma_r_q=sigma_6*q**3+sigma_7*q**2+sigma_8*q+sigma_9
    mu_r=(mu_r_lambda_s+mu_r_q)/2.
    sigma_r=(sigma_r_lambda_s**2+sigma_r_q**2)**(0.5)
    lambda_a_marginalized=lambda_a+np.random.normal(loc=mu_r,scale=sigma_r)

    return lambda_a_marginalized

def convert_lambda_s_to_lambda_a_numpy(lambda_s,q):
    
    #define bmatrix and cmatrix
    bmatrix=np.array([[-14.4,14.45],[31.36,-32.25],[-22.44,20.35]])
    cmatrix=np.array([[-15.25,15.37],[37.33,-43.20],[-29.93,35.18]])
    
    #define other quantities for computing the binary love relation
    n=0.743
    alpha=1
    Fn=(1-q**(10/(3-n)))/(1+q**(10/(3-n)))
    q_arr = q**np.arange(1,3)
    lambda_arr = lambda_s**(-np.arange(1,4)/5.)
    #compute the sum in the numerator and denominator of binary love relation
    num_sum= np.dot(lambda_arr,bmatrix.dot(q_arr))
    den_sum= np.dot(lambda_arr,cmatrix.dot(q_arr))
    #compute binary love relation
    lambda_a=Fn*(1+num_sum)/(1+den_sum)*lambda_s**(alpha)
    
    return lambda_a


def convert_lambda_s_to_lambda_a_marginalized_numpy(lambda_s,q):
    """
    Marginalized Binary Love relations. 
    
    Lambda_s = (Lambda_1 + Lambda_2)/2
    Lambda_a = (Lambda_1 - Lambda_2)/2
   
    The Binary Love relations give us Lambda_a(Lambda_s).

    This relation has some error \delta. 
    Following [1], we marginalize over this error by sampling from 
    a Gaussian distributions with mean/variance
    equal to the mean/variance of the distribution of predictions from the 
    binary Love relation for different EOS.
    Thus this function will give you a different answer if you call it multiple
    times with the same input, unless you manually reset the random seed every
    time you call this function.

    [1] https://arxiv.org/abs/1903.03909

    DOI: https://doi.org/10.1103/PhysRevD.99.083016
    """ 
    
    lambda_a = convert_lambda_s_to_lambda_a_numpy(lambda_s,q)

    #define coefficients mu's and sigma's
    mu_1=3.509e-3
    mu_2=9.351e-1
    mu_3=-18.07
    mu_4=27.56
    mu_5=-10.10
    sigma_1=-2.074e-7
    sigma_2=1.492e-3
    sigma_3=-4.891e-2
    sigma_4=8.207e-1
    sigma_5=-1.308
    sigma_6=-63.76
    sigma_7=11.14
    sigma_8=75.25
    sigma_9=-23.69 

    #construct the marginalized mean and standard deviation 
    mu_r_lambda_s=mu_1*lambda_s+mu_2
    mu_r_q=mu_3*q**2+mu_4*q+mu_5
    sigma_r_lambda_s=(sigma_1*lambda_s**(5./2)+sigma_2*lambda_s**(3./2) +sigma_3*lambda_s+sigma_4*lambda_s**(1./2)+sigma_5)
    sigma_r_q=sigma_6*q**3+sigma_7*q**2+sigma_8*q+sigma_9
    mu_r=(mu_r_lambda_s+mu_r_q)/2.
    sigma_r=(sigma_r_lambda_s**2+sigma_r_q**2)**(0.5)
    #print(lambda_a)
    lambda_a_marginalized=lambda_a+np.random.normal(loc=mu_r,scale=sigma_r)
    
    return lambda_a_marginalized
