#!/usr/bin/env python
# coding: utf-8

# # Functions

# ## Vanilla Euro

# In[1]:


# =============================================================================
# This file is especially for visualization of option price/Greeks w.r.t its parameters
# Mainly based on Quantlib
# =============================================================================
import QuantLib as ql

def vanilla_euro_call(S,K,r,q,tau,sigma,**kwargs):
    # specify timeframe
    tau_day=int(round(tau*365))
    start=ql.Date(1,1,2000)
    ql.Settings.instance().evaluationDate=start
    maturity=start+tau_day
    day_count = ql.Actual365Fixed()
    calendar = ql.UnitedStates()
    # specify option type
    option_type = ql.Option.Call
    payoff = ql.PlainVanillaPayoff(option_type, K)
    exercise = ql.EuropeanExercise(maturity)
    option=ql.VanillaOption(payoff, exercise)
    # specify pricing method
    spot_handle = ql.QuoteHandle(
    ql.SimpleQuote(S)
    )
    flat_ts = ql.YieldTermStructureHandle(
        ql.FlatForward(start, r, day_count)
    )
    dividend_yield = ql.YieldTermStructureHandle(
        ql.FlatForward(start, q, day_count)
    )
    flat_vol_ts = ql.BlackVolTermStructureHandle(
        ql.BlackConstantVol(start, calendar, sigma, day_count)
    )
    bsm_process = ql.BlackScholesMertonProcess(spot_handle, 
                                           dividend_yield, 
                                           flat_ts, 
                                           flat_vol_ts)
    option.setPricingEngine(ql.AnalyticEuropeanEngine(bsm_process))
    #setup is finished
    # return option.NPV(),option.delta(),option.gamma(),option.vega(),option.theta(),option.rho()
    return {
        'PV':option.NPV(),
        'delta':option.delta(),
        'gamma':option.gamma(),
        'vega':option.vega(),
        'theta':option.theta(),
        'rho':option.rho(),
        }

def vanilla_euro_put(S,K,r,q,tau,sigma,**kwargs):
    # specify timeframe
    tau_day=int(round(tau*365))
    start=ql.Date(1,1,2000)
    ql.Settings.instance().evaluationDate=start
    maturity=start+tau_day
    day_count = ql.Actual365Fixed()
    calendar = ql.UnitedStates()
    # specify option type
    option_type = ql.Option.Put
    payoff = ql.PlainVanillaPayoff(option_type, K)
    exercise = ql.EuropeanExercise(maturity)
    option=ql.VanillaOption(payoff, exercise)
    # specify pricing method
    spot_handle = ql.QuoteHandle(
    ql.SimpleQuote(S)
    )
    flat_ts = ql.YieldTermStructureHandle(
        ql.FlatForward(start, r, day_count)
    )
    dividend_yield = ql.YieldTermStructureHandle(
        ql.FlatForward(start, q, day_count)
    )
    flat_vol_ts = ql.BlackVolTermStructureHandle(
        ql.BlackConstantVol(start, calendar, sigma, day_count)
    )
    bsm_process = ql.BlackScholesMertonProcess(spot_handle, 
                                           dividend_yield, 
                                           flat_ts, 
                                           flat_vol_ts)
    option.setPricingEngine(ql.AnalyticEuropeanEngine(bsm_process))
    #setup is finished
    # return option.NPV(),option.delta(),option.gamma(),option.vega(),option.theta(),option.rho()
    return {
        'PV':option.NPV(),
        'delta':option.delta(),
        'gamma':option.gamma(),
        'vega':option.vega(),
        'theta':option.theta(),
        'rho':option.rho(),
        }


# ## Vanilla Ame

# In[2]:


def _vanilla_ame_call(S,K,r,q,tau,sigma):
    # specify timeframe
    tau_day=int(round(tau*365))
    start=ql.Date(1,1,2000)
    ql.Settings.instance().evaluationDate=start
    maturity=start+tau_day
    day_count = ql.Actual365Fixed()
    calendar = ql.UnitedStates()
    # specify option type
    option_type = ql.Option.Call
    payoff = ql.PlainVanillaPayoff(option_type, K)
    exercise = ql.AmericanExercise(start,maturity)
    option=ql.VanillaOption(payoff, exercise)
    # specify pricing method
    spot_handle = ql.QuoteHandle(
    ql.SimpleQuote(S)
    )
    flat_ts = ql.YieldTermStructureHandle(
        ql.FlatForward(start, r, day_count)
    )
    dividend_yield = ql.YieldTermStructureHandle(
        ql.FlatForward(start, q, day_count)
    )
    flat_vol_ts = ql.BlackVolTermStructureHandle(
        ql.BlackConstantVol(start, calendar, sigma, day_count)
    )
    bsm_process = ql.BlackScholesMertonProcess(spot_handle, 
                                           dividend_yield, 
                                           flat_ts, 
                                           flat_vol_ts)
    steps=tau_day*2
    option.setPricingEngine(ql.BinomialVanillaEngine(bsm_process,'crr',steps))
    #setup is finished
    return option

def vanilla_ame_call(S,K,r,q,tau,sigma,**kwargs):
    option=_vanilla_ame_call(S,K,r,q,tau,sigma)
    epsilon=1e-8
    option_sigam_upper=_vanilla_ame_call(S,K,r,q,tau,sigma+epsilon)
    option_sigam_lower=_vanilla_ame_call(S,K,r,q,tau,sigma-epsilon)
    option_r_upper=_vanilla_ame_call(S,K,r+epsilon,q,tau,sigma)
    option_r_lower=_vanilla_ame_call(S,K,r-epsilon,q,tau,sigma)
    return {
        'PV':option.NPV(),
        'delta':option.delta(),
        'gamma':option.gamma(),
        'vega':(option_sigam_upper.NPV()-option_sigam_lower.NPV())/(2*epsilon),
        'theta':option.theta(),
        'rho':(option_r_upper.NPV()-option_r_lower.NPV())/(2*epsilon),
        }

def _vanilla_ame_put(S,K,r,q,tau,sigma):
    # specify timeframe
    tau_day=int(round(tau*365))
    start=ql.Date(1,1,2000)
    ql.Settings.instance().evaluationDate=start
    maturity=start+tau_day
    day_count = ql.Actual365Fixed()
    calendar = ql.UnitedStates()
    # specify option type
    option_type = ql.Option.Put
    payoff = ql.PlainVanillaPayoff(option_type, K)
    exercise = ql.AmericanExercise(start,maturity)
    option=ql.VanillaOption(payoff, exercise)
    # specify pricing method
    spot_handle = ql.QuoteHandle(
    ql.SimpleQuote(S)
    )
    flat_ts = ql.YieldTermStructureHandle(
        ql.FlatForward(start, r, day_count)
    )
    dividend_yield = ql.YieldTermStructureHandle(
        ql.FlatForward(start, q, day_count)
    )
    flat_vol_ts = ql.BlackVolTermStructureHandle(
        ql.BlackConstantVol(start, calendar, sigma, day_count)
    )
    bsm_process = ql.BlackScholesMertonProcess(spot_handle, 
                                           dividend_yield, 
                                           flat_ts, 
                                           flat_vol_ts)
    steps=tau_day*2
    option.setPricingEngine(ql.BinomialVanillaEngine(bsm_process,'crr',steps))
    #setup is finished
    return option
def vanilla_ame_put(S,K,r,q,tau,sigma,**kwargs):
    option=_vanilla_ame_put(S,K,r,q,tau,sigma)
    epsilon=1e-8
    option_sigam_upper=_vanilla_ame_put(S,K,r,q,tau,sigma+epsilon)
    option_sigam_lower=_vanilla_ame_put(S,K,r,q,tau,sigma-epsilon)
    option_r_upper=_vanilla_ame_put(S,K,r+epsilon,q,tau,sigma)
    option_r_lower=_vanilla_ame_put(S,K,r-epsilon,q,tau,sigma)
    return {
        'PV':option.NPV(),
        'delta':option.delta(),
        'gamma':option.gamma(),
        'vega':(option_sigam_upper.NPV()-option_sigam_lower.NPV())/(2*epsilon),
        'theta':option.theta(),
        'rho':(option_r_upper.NPV()-option_r_lower.NPV())/(2*epsilon),
        }


# ## Barrier

# In[26]:


def _barrier_upin_call(S,K,X,r,q,tau,sigma):
    # specify timeframe
    tau_day=int(round(tau*365))
    start=ql.Date(1,1,2000)
    ql.Settings.instance().evaluationDate=start
    maturity=start+tau_day
    day_count = ql.Actual365Fixed()
    calendar = ql.UnitedStates()
    # specify option type
    option_type = ql.Option.Call
    payoff = ql.PlainVanillaPayoff(option_type, K)
    exercise = ql.EuropeanExercise(maturity)
    option=ql.BarrierOption(ql.Barrier.UpIn, X, 0, payoff, exercise)
    # specify pricing method
    spot_handle = ql.QuoteHandle(
    ql.SimpleQuote(S)
    )
    flat_ts = ql.YieldTermStructureHandle(
        ql.FlatForward(start, r, day_count)
    )
    dividend_yield = ql.YieldTermStructureHandle(
        ql.FlatForward(start, q, day_count)
    )
    flat_vol_ts = ql.BlackVolTermStructureHandle(
        ql.BlackConstantVol(start, calendar, sigma, day_count)
    )
    bsm_process = ql.BlackScholesMertonProcess(spot_handle, 
                                           dividend_yield, 
                                           flat_ts, 
                                           flat_vol_ts)
    option.setPricingEngine(ql.AnalyticBarrierEngine(bsm_process))
    return option

def barrier_upin_call(S,K,X,r,q,tau,sigma,**kwargs):
    option=_barrier_upin_call(S,K,X,r,q,tau,sigma)
    epsilon=1e-4
    option_S_upper=_barrier_upin_call(S+epsilon,K,X,r,q,tau,sigma)
    option_S_lower=_barrier_upin_call(S-epsilon,K,X,r,q,tau,sigma)
    option_sigam_upper=_barrier_upin_call(S,K,X,r,q,tau,sigma+epsilon)
    option_sigam_lower=_barrier_upin_call(S,K,X,r,q,tau,sigma-epsilon)
    option_tau_upper=_barrier_upin_call(S,K,X,r,q,tau+100*epsilon,sigma)
    option_tau_lower=_barrier_upin_call(S,K,X,r,q,tau-100*epsilon,sigma)
    option_r_upper=_barrier_upin_call(S,K,X,r+epsilon,q,tau,sigma)
    option_r_lower=_barrier_upin_call(S,K,X,r-epsilon,q,tau,sigma)
    return {
        'PV':option.NPV(),
        'delta':(option_S_upper.NPV()-option_S_lower.NPV())/(2*epsilon),
        'gamma':(option_S_upper.NPV()-2*option.NPV()+option_S_lower.NPV())/(epsilon**2),
        'vega':(option_sigam_upper.NPV()-option_sigam_lower.NPV())/(2*epsilon),
        'theta':(option_tau_upper.NPV()-option_tau_lower.NPV())/(2*100*epsilon),
        'rho':(option_r_upper.NPV()-option_r_lower.NPV())/(2*epsilon),
        }

def _barrier_upin_put(S,K,X,r,q,tau,sigma):
    # specify timeframe
    tau_day=int(round(tau*365))
    start=ql.Date(1,1,2000)
    ql.Settings.instance().evaluationDate=start
    maturity=start+tau_day
    day_count = ql.Actual365Fixed()
    calendar = ql.UnitedStates()
    # specify option type
    option_type = ql.Option.Put
    payoff = ql.PlainVanillaPayoff(option_type, K)
    exercise = ql.EuropeanExercise(maturity)
    option=ql.BarrierOption(ql.Barrier.UpIn, X, 0, payoff, exercise)
    # specify pricing method
    spot_handle = ql.QuoteHandle(
    ql.SimpleQuote(S)
    )
    flat_ts = ql.YieldTermStructureHandle(
        ql.FlatForward(start, r, day_count)
    )
    dividend_yield = ql.YieldTermStructureHandle(
        ql.FlatForward(start, q, day_count)
    )
    flat_vol_ts = ql.BlackVolTermStructureHandle(
        ql.BlackConstantVol(start, calendar, sigma, day_count)
    )
    bsm_process = ql.BlackScholesMertonProcess(spot_handle, 
                                           dividend_yield, 
                                           flat_ts, 
                                           flat_vol_ts)
    option.setPricingEngine(ql.AnalyticBarrierEngine(bsm_process))
    return option

def barrier_upin_put(S,K,X,r,q,tau,sigma,**kwargs):
    option=_barrier_upin_put(S,K,X,r,q,tau,sigma)
    epsilon=1e-4
    option_S_upper=_barrier_upin_put(S+epsilon,K,X,r,q,tau,sigma)
    option_S_lower=_barrier_upin_put(S-epsilon,K,X,r,q,tau,sigma)
    option_sigam_upper=_barrier_upin_put(S,K,X,r,q,tau,sigma+epsilon)
    option_sigam_lower=_barrier_upin_put(S,K,X,r,q,tau,sigma-epsilon)
    option_tau_upper=_barrier_upin_put(S,K,X,r,q,tau+100*epsilon,sigma)
    option_tau_lower=_barrier_upin_put(S,K,X,r,q,tau-100*epsilon,sigma)
    option_r_upper=_barrier_upin_put(S,K,X,r+epsilon,q,tau,sigma)
    option_r_lower=_barrier_upin_put(S,K,X,r-epsilon,q,tau,sigma)
    return {
        'PV':option.NPV(),
        'delta':(option_S_upper.NPV()-option_S_lower.NPV())/(2*epsilon),
        'gamma':(option_S_upper.NPV()-2*option.NPV()+option_S_lower.NPV())/(epsilon**2),
        'vega':(option_sigam_upper.NPV()-option_sigam_lower.NPV())/(2*epsilon),
        'theta':(option_tau_upper.NPV()-option_tau_lower.NPV())/(2*100*epsilon),
        'rho':(option_r_upper.NPV()-option_r_lower.NPV())/(2*epsilon),
        }

def _barrier_upout_call(S,K,X,r,q,tau,sigma):
    # specify timeframe
    tau_day=int(round(tau*365))
    start=ql.Date(1,1,2000)
    ql.Settings.instance().evaluationDate=start
    maturity=start+tau_day
    day_count = ql.Actual365Fixed()
    calendar = ql.UnitedStates()
    # specify option type
    option_type = ql.Option.Call
    payoff = ql.PlainVanillaPayoff(option_type, K)
    exercise = ql.EuropeanExercise(maturity)
    option=ql.BarrierOption(ql.Barrier.UpOut, X, 0, payoff, exercise)
    # specify pricing method
    spot_handle = ql.QuoteHandle(
    ql.SimpleQuote(S)
    )
    flat_ts = ql.YieldTermStructureHandle(
        ql.FlatForward(start, r, day_count)
    )
    dividend_yield = ql.YieldTermStructureHandle(
        ql.FlatForward(start, q, day_count)
    )
    flat_vol_ts = ql.BlackVolTermStructureHandle(
        ql.BlackConstantVol(start, calendar, sigma, day_count)
    )
    bsm_process = ql.BlackScholesMertonProcess(spot_handle, 
                                           dividend_yield, 
                                           flat_ts, 
                                           flat_vol_ts)
    option.setPricingEngine(ql.AnalyticBarrierEngine(bsm_process))
    return option

def barrier_upout_call(S,K,X,r,q,tau,sigma,**kwargs):
    option=_barrier_upout_call(S,K,X,r,q,tau,sigma)
    epsilon=1e-4
    option_S_upper=_barrier_upout_call(S+epsilon,K,X,r,q,tau,sigma)
    option_S_lower=_barrier_upout_call(S-epsilon,K,X,r,q,tau,sigma)
    option_sigam_upper=_barrier_upout_call(S,K,X,r,q,tau,sigma+epsilon)
    option_sigam_lower=_barrier_upout_call(S,K,X,r,q,tau,sigma-epsilon)
    option_tau_upper=_barrier_upout_call(S,K,X,r,q,tau+100*epsilon,sigma)
    option_tau_lower=_barrier_upout_call(S,K,X,r,q,tau-100*epsilon,sigma)
    option_r_upper=_barrier_upout_call(S,K,X,r+epsilon,q,tau,sigma)
    option_r_lower=_barrier_upout_call(S,K,X,r-epsilon,q,tau,sigma)
    return {
        'PV':option.NPV(),
        'delta':(option_S_upper.NPV()-option_S_lower.NPV())/(2*epsilon),
        'gamma':(option_S_upper.NPV()-2*option.NPV()+option_S_lower.NPV())/(epsilon**2),
        'vega':(option_sigam_upper.NPV()-option_sigam_lower.NPV())/(2*epsilon),
        'theta':(option_tau_upper.NPV()-option_tau_lower.NPV())/(2*100*epsilon),
        'rho':(option_r_upper.NPV()-option_r_lower.NPV())/(2*epsilon),
        }

def _barrier_upout_put(S,K,X,r,q,tau,sigma):
    # specify timeframe
    tau_day=int(round(tau*365))
    start=ql.Date(1,1,2000)
    ql.Settings.instance().evaluationDate=start
    maturity=start+tau_day
    day_count = ql.Actual365Fixed()
    calendar = ql.UnitedStates()
    # specify option type
    option_type = ql.Option.Put
    payoff = ql.PlainVanillaPayoff(option_type, K)
    exercise = ql.EuropeanExercise(maturity)
    option=ql.BarrierOption(ql.Barrier.UpOut, X, 0, payoff, exercise)
    # specify pricing method
    spot_handle = ql.QuoteHandle(
    ql.SimpleQuote(S)
    )
    flat_ts = ql.YieldTermStructureHandle(
        ql.FlatForward(start, r, day_count)
    )
    dividend_yield = ql.YieldTermStructureHandle(
        ql.FlatForward(start, q, day_count)
    )
    flat_vol_ts = ql.BlackVolTermStructureHandle(
        ql.BlackConstantVol(start, calendar, sigma, day_count)
    )
    bsm_process = ql.BlackScholesMertonProcess(spot_handle, 
                                           dividend_yield, 
                                           flat_ts, 
                                           flat_vol_ts)
    option.setPricingEngine(ql.AnalyticBarrierEngine(bsm_process))
    return option

def barrier_upout_put(S,K,X,r,q,tau,sigma,**kwargs):
    option=_barrier_upout_put(S,K,X,r,q,tau,sigma)
    epsilon=1e-4
    option_S_upper=_barrier_upout_put(S+epsilon,K,X,r,q,tau,sigma)
    option_S_lower=_barrier_upout_put(S-epsilon,K,X,r,q,tau,sigma)
    option_sigam_upper=_barrier_upout_put(S,K,X,r,q,tau,sigma+epsilon)
    option_sigam_lower=_barrier_upout_put(S,K,X,r,q,tau,sigma-epsilon)
    option_tau_upper=_barrier_upout_put(S,K,X,r,q,tau+100*epsilon,sigma)
    option_tau_lower=_barrier_upout_put(S,K,X,r,q,tau-100*epsilon,sigma)
    option_r_upper=_barrier_upout_put(S,K,X,r+epsilon,q,tau,sigma)
    option_r_lower=_barrier_upout_put(S,K,X,r-epsilon,q,tau,sigma)
    return {
        'PV':option.NPV(),
        'delta':(option_S_upper.NPV()-option_S_lower.NPV())/(2*epsilon),
        'gamma':(option_S_upper.NPV()-2*option.NPV()+option_S_lower.NPV())/(epsilon**2),
        'vega':(option_sigam_upper.NPV()-option_sigam_lower.NPV())/(2*epsilon),
        'theta':(option_tau_upper.NPV()-option_tau_lower.NPV())/(2*100*epsilon),
        'rho':(option_r_upper.NPV()-option_r_lower.NPV())/(2*epsilon),
        }

def _barrier_downin_call(S,K,X,r,q,tau,sigma):
    # specify timeframe
    tau_day=int(round(tau*365))
    start=ql.Date(1,1,2000)
    ql.Settings.instance().evaluationDate=start
    maturity=start+tau_day
    day_count = ql.Actual365Fixed()
    calendar = ql.UnitedStates()
    # specify option type
    option_type = ql.Option.Call
    payoff = ql.PlainVanillaPayoff(option_type, K)
    exercise = ql.EuropeanExercise(maturity)
    option=ql.BarrierOption(ql.Barrier.DownIn, X, 0, payoff, exercise)
    # specify pricing method
    spot_handle = ql.QuoteHandle(
    ql.SimpleQuote(S)
    )
    flat_ts = ql.YieldTermStructureHandle(
        ql.FlatForward(start, r, day_count)
    )
    dividend_yield = ql.YieldTermStructureHandle(
        ql.FlatForward(start, q, day_count)
    )
    flat_vol_ts = ql.BlackVolTermStructureHandle(
        ql.BlackConstantVol(start, calendar, sigma, day_count)
    )
    bsm_process = ql.BlackScholesMertonProcess(spot_handle, 
                                           dividend_yield, 
                                           flat_ts, 
                                           flat_vol_ts)
    option.setPricingEngine(ql.AnalyticBarrierEngine(bsm_process))
    return option

def barrier_downin_call(S,K,X,r,q,tau,sigma,**kwargs):
    option=_barrier_downin_call(S,K,X,r,q,tau,sigma)
    epsilon=1e-4
    option_S_upper=_barrier_downin_call(S+epsilon,K,X,r,q,tau,sigma)
    option_S_lower=_barrier_downin_call(S-epsilon,K,X,r,q,tau,sigma)
    option_sigam_upper=_barrier_downin_call(S,K,X,r,q,tau,sigma+epsilon)
    option_sigam_lower=_barrier_downin_call(S,K,X,r,q,tau,sigma-epsilon)
    option_tau_upper=_barrier_downin_call(S,K,X,r,q,tau+100*epsilon,sigma)
    option_tau_lower=_barrier_downin_call(S,K,X,r,q,tau-100*epsilon,sigma)
    option_r_upper=_barrier_downin_call(S,K,X,r+epsilon,q,tau,sigma)
    option_r_lower=_barrier_downin_call(S,K,X,r-epsilon,q,tau,sigma)
    return {
        'PV':option.NPV(),
        'delta':(option_S_upper.NPV()-option_S_lower.NPV())/(2*epsilon),
        'gamma':(option_S_upper.NPV()-2*option.NPV()+option_S_lower.NPV())/(epsilon**2),
        'vega':(option_sigam_upper.NPV()-option_sigam_lower.NPV())/(2*epsilon),
        'theta':(option_tau_upper.NPV()-option_tau_lower.NPV())/(2*100*epsilon),
        'rho':(option_r_upper.NPV()-option_r_lower.NPV())/(2*epsilon),
        }

def _barrier_downin_put(S,K,X,r,q,tau,sigma):
    # specify timeframe
    tau_day=int(round(tau*365))
    start=ql.Date(1,1,2000)
    ql.Settings.instance().evaluationDate=start
    maturity=start+tau_day
    day_count = ql.Actual365Fixed()
    calendar = ql.UnitedStates()
    # specify option type
    option_type = ql.Option.Put
    payoff = ql.PlainVanillaPayoff(option_type, K)
    exercise = ql.EuropeanExercise(maturity)
    option=ql.BarrierOption(ql.Barrier.DownIn, X, 0, payoff, exercise)
    # specify pricing method
    spot_handle = ql.QuoteHandle(
    ql.SimpleQuote(S)
    )
    flat_ts = ql.YieldTermStructureHandle(
        ql.FlatForward(start, r, day_count)
    )
    dividend_yield = ql.YieldTermStructureHandle(
        ql.FlatForward(start, q, day_count)
    )
    flat_vol_ts = ql.BlackVolTermStructureHandle(
        ql.BlackConstantVol(start, calendar, sigma, day_count)
    )
    bsm_process = ql.BlackScholesMertonProcess(spot_handle, 
                                           dividend_yield, 
                                           flat_ts, 
                                           flat_vol_ts)
    option.setPricingEngine(ql.AnalyticBarrierEngine(bsm_process))
    return option

def barrier_downin_put(S,K,X,r,q,tau,sigma,**kwargs):
    option=_barrier_downin_put(S,K,X,r,q,tau,sigma)
    epsilon=1e-4
    option_S_upper=_barrier_downin_put(S+epsilon,K,X,r,q,tau,sigma)
    option_S_lower=_barrier_downin_put(S-epsilon,K,X,r,q,tau,sigma)
    option_sigam_upper=_barrier_downin_put(S,K,X,r,q,tau,sigma+epsilon)
    option_sigam_lower=_barrier_downin_put(S,K,X,r,q,tau,sigma-epsilon)
    option_tau_upper=_barrier_downin_put(S,K,X,r,q,tau+100*epsilon,sigma)
    option_tau_lower=_barrier_downin_put(S,K,X,r,q,tau-100*epsilon,sigma)
    option_r_upper=_barrier_downin_put(S,K,X,r+epsilon,q,tau,sigma)
    option_r_lower=_barrier_downin_put(S,K,X,r-epsilon,q,tau,sigma)
    return {
        'PV':option.NPV(),
        'delta':(option_S_upper.NPV()-option_S_lower.NPV())/(2*epsilon),
        'gamma':(option_S_upper.NPV()-2*option.NPV()+option_S_lower.NPV())/(epsilon**2),
        'vega':(option_sigam_upper.NPV()-option_sigam_lower.NPV())/(2*epsilon),
        'theta':(option_tau_upper.NPV()-option_tau_lower.NPV())/(2*100*epsilon),
        'rho':(option_r_upper.NPV()-option_r_lower.NPV())/(2*epsilon),
        }

def _barrier_downout_call(S,K,X,r,q,tau,sigma):
    # specify timeframe
    tau_day=int(round(tau*365))
    start=ql.Date(1,1,2000)
    ql.Settings.instance().evaluationDate=start
    maturity=start+tau_day
    day_count = ql.Actual365Fixed()
    calendar = ql.UnitedStates()
    # specify option type
    option_type = ql.Option.Call
    payoff = ql.PlainVanillaPayoff(option_type, K)
    exercise = ql.EuropeanExercise(maturity)
    option=ql.BarrierOption(ql.Barrier.DownOut, X, 0, payoff, exercise)
    # specify pricing method
    spot_handle = ql.QuoteHandle(
    ql.SimpleQuote(S)
    )
    flat_ts = ql.YieldTermStructureHandle(
        ql.FlatForward(start, r, day_count)
    )
    dividend_yield = ql.YieldTermStructureHandle(
        ql.FlatForward(start, q, day_count)
    )
    flat_vol_ts = ql.BlackVolTermStructureHandle(
        ql.BlackConstantVol(start, calendar, sigma, day_count)
    )
    bsm_process = ql.BlackScholesMertonProcess(spot_handle, 
                                           dividend_yield, 
                                           flat_ts, 
                                           flat_vol_ts)
    option.setPricingEngine(ql.AnalyticBarrierEngine(bsm_process))
    return option

def barrier_downout_call(S,K,X,r,q,tau,sigma,**kwargs):
    option=_barrier_downout_call(S,K,X,r,q,tau,sigma)
    epsilon=1e-4
    option_S_upper=_barrier_downout_call(S+epsilon,K,X,r,q,tau,sigma)
    option_S_lower=_barrier_downout_call(S-epsilon,K,X,r,q,tau,sigma)
    option_sigam_upper=_barrier_downout_call(S,K,X,r,q,tau,sigma+epsilon)
    option_sigam_lower=_barrier_downout_call(S,K,X,r,q,tau,sigma-epsilon)
    option_tau_upper=_barrier_downout_call(S,K,X,r,q,tau+100*epsilon,sigma)
    option_tau_lower=_barrier_downout_call(S,K,X,r,q,tau-100*epsilon,sigma)
    option_r_upper=_barrier_downout_call(S,K,X,r+epsilon,q,tau,sigma)
    option_r_lower=_barrier_downout_call(S,K,X,r-epsilon,q,tau,sigma)
    return {
        'PV':option.NPV(),
        'delta':(option_S_upper.NPV()-option_S_lower.NPV())/(2*epsilon),
        'gamma':(option_S_upper.NPV()-2*option.NPV()+option_S_lower.NPV())/(epsilon**2),
        'vega':(option_sigam_upper.NPV()-option_sigam_lower.NPV())/(2*epsilon),
        'theta':(option_tau_upper.NPV()-option_tau_lower.NPV())/(2*100*epsilon),
        'rho':(option_r_upper.NPV()-option_r_lower.NPV())/(2*epsilon),
        }

def _barrier_downout_put(S,K,X,r,q,tau,sigma):
    # specify timeframe
    tau_day=int(round(tau*365))
    start=ql.Date(1,1,2000)
    ql.Settings.instance().evaluationDate=start
    maturity=start+tau_day
    day_count = ql.Actual365Fixed()
    calendar = ql.UnitedStates()
    # specify option type
    option_type = ql.Option.Put
    payoff = ql.PlainVanillaPayoff(option_type, K)
    exercise = ql.EuropeanExercise(maturity)
    option=ql.BarrierOption(ql.Barrier.DownOut, X, 0, payoff, exercise)
    # specify pricing method
    spot_handle = ql.QuoteHandle(
    ql.SimpleQuote(S)
    )
    flat_ts = ql.YieldTermStructureHandle(
        ql.FlatForward(start, r, day_count)
    )
    dividend_yield = ql.YieldTermStructureHandle(
        ql.FlatForward(start, q, day_count)
    )
    flat_vol_ts = ql.BlackVolTermStructureHandle(
        ql.BlackConstantVol(start, calendar, sigma, day_count)
    )
    bsm_process = ql.BlackScholesMertonProcess(spot_handle, 
                                           dividend_yield, 
                                           flat_ts, 
                                           flat_vol_ts)
    option.setPricingEngine(ql.AnalyticBarrierEngine(bsm_process))
    return option

def barrier_downout_put(S,K,X,r,q,tau,sigma,**kwargs):
    option=_barrier_downout_put(S,K,X,r,q,tau,sigma)
    epsilon=1e-4
    option_S_upper=_barrier_downout_put(S+epsilon,K,X,r,q,tau,sigma)
    option_S_lower=_barrier_downout_put(S-epsilon,K,X,r,q,tau,sigma)
    option_sigam_upper=_barrier_downout_put(S,K,X,r,q,tau,sigma+epsilon)
    option_sigam_lower=_barrier_downout_put(S,K,X,r,q,tau,sigma-epsilon)
    option_tau_upper=_barrier_downout_put(S,K,X,r,q,tau+100*epsilon,sigma)
    option_tau_lower=_barrier_downout_put(S,K,X,r,q,tau-100*epsilon,sigma)
    option_r_upper=_barrier_downout_put(S,K,X,r+epsilon,q,tau,sigma)
    option_r_lower=_barrier_downout_put(S,K,X,r-epsilon,q,tau,sigma)
    return {
        'PV':option.NPV(),
        'delta':(option_S_upper.NPV()-option_S_lower.NPV())/(2*epsilon),
        'gamma':(option_S_upper.NPV()-2*option.NPV()+option_S_lower.NPV())/(epsilon**2),
        'vega':(option_sigam_upper.NPV()-option_sigam_lower.NPV())/(2*epsilon),
        'theta':(option_tau_upper.NPV()-option_tau_lower.NPV())/(2*100*epsilon),
        'rho':(option_r_upper.NPV()-option_r_lower.NPV())/(2*epsilon),
        }


# ## Asian

# In[23]:


def asian_geometric_call(S,K,r,q,tau,sigma,freq='1M',**kwargs):
    # specify timeframe
    tau_day=int(round(tau*365))
    start=ql.Date(1,1,2000)
    ql.Settings.instance().evaluationDate=start
    maturity=start+tau_day
    day_count = ql.Actual365Fixed()
    calendar = ql.UnitedStates()
    period=ql.Period(freq)
    asianFutureFixingDates = [start+period]
    while asianFutureFixingDates[-1]<=maturity:
        asianFutureFixingDates.append(asianFutureFixingDates[-1]+period)
    asianFutureFixingDates.pop(-1)
    # specify option type
    option_type = ql.Option.Call
    payoff = ql.PlainVanillaPayoff(option_type, K)
    exercise = ql.EuropeanExercise(maturity)
    option = ql.DiscreteAveragingAsianOption(ql.Average().Geometric,
                                             1, 0, 
                                             asianFutureFixingDates, 
                                             payoff, 
                                             exercise)
    # specify pricing method
    spot_handle = ql.QuoteHandle(
    ql.SimpleQuote(S)
    )
    flat_ts = ql.YieldTermStructureHandle(
        ql.FlatForward(start, r, day_count)
    )
    dividend_yield = ql.YieldTermStructureHandle(
        ql.FlatForward(start, q, day_count)
    )
    flat_vol_ts = ql.BlackVolTermStructureHandle(
        ql.BlackConstantVol(start, calendar, sigma, day_count)
    )
    bsm_process = ql.BlackScholesMertonProcess(spot_handle, 
                                           dividend_yield, 
                                           flat_ts, 
                                           flat_vol_ts)
    option.setPricingEngine(ql.AnalyticDiscreteGeometricAveragePriceAsianEngine(bsm_process))
    return {
    'PV':option.NPV(),
    'delta':option.delta(),
    'gamma':option.gamma(),
    'vega':option.vega(),
    'theta':option.theta(),
    'rho':option.rho(),
    }

def asian_geometric_put(S,K,r,q,tau,sigma,freq='1M',**kwargs):
    # specify timeframe
    tau_day=int(round(tau*365))
    start=ql.Date(1,1,2000)
    ql.Settings.instance().evaluationDate=start
    maturity=start+tau_day
    day_count = ql.Actual365Fixed()
    calendar = ql.UnitedStates()
    period=ql.Period(freq)
    asianFutureFixingDates = [start+period]
    while asianFutureFixingDates[-1]<=maturity:
        asianFutureFixingDates.append(asianFutureFixingDates[-1]+period)
    asianFutureFixingDates.pop(-1)
    # specify option type
    option_type = ql.Option.Put
    payoff = ql.PlainVanillaPayoff(option_type, K)
    exercise = ql.EuropeanExercise(maturity)
    option = ql.DiscreteAveragingAsianOption(ql.Average().Geometric,
                                             1, 0, 
                                             asianFutureFixingDates, 
                                             payoff, 
                                             exercise)
    # specify pricing method
    spot_handle = ql.QuoteHandle(
    ql.SimpleQuote(S)
    )
    flat_ts = ql.YieldTermStructureHandle(
        ql.FlatForward(start, r, day_count)
    )
    dividend_yield = ql.YieldTermStructureHandle(
        ql.FlatForward(start, q, day_count)
    )
    flat_vol_ts = ql.BlackVolTermStructureHandle(
        ql.BlackConstantVol(start, calendar, sigma, day_count)
    )
    bsm_process = ql.BlackScholesMertonProcess(spot_handle, 
                                           dividend_yield, 
                                           flat_ts, 
                                           flat_vol_ts)
    option.setPricingEngine(ql.AnalyticDiscreteGeometricAveragePriceAsianEngine(bsm_process))
    return {
    'PV':option.NPV(),
    'delta':option.delta(),
    'gamma':option.gamma(),
    'vega':option.vega(),
    'theta':option.theta(),
    'rho':option.rho(),
    }

def _asian_arithmetic_call(S,K,r,q,tau,sigma,freq='1M'):
    # specify timeframe
    tau_day=int(round(tau*365))
    start=ql.Date(1,1,2000)
    ql.Settings.instance().evaluationDate=start
    maturity=start+tau_day
    day_count = ql.Actual365Fixed()
    calendar = ql.UnitedStates()
    period=ql.Period(freq)
    asianFutureFixingDates = [start+period]
    while asianFutureFixingDates[-1]<=maturity:
        asianFutureFixingDates.append(asianFutureFixingDates[-1]+period)
    asianFutureFixingDates.pop(-1)
    # specify option type
    option_type = ql.Option.Call
    payoff = ql.PlainVanillaPayoff(option_type, K)
    exercise = ql.EuropeanExercise(maturity)
    option = ql.DiscreteAveragingAsianOption(ql.Average().Arithmetic,
                                             0, 0, 
                                             asianFutureFixingDates, 
                                             payoff, 
                                             exercise)
    # specify pricing method
    spot_handle = ql.QuoteHandle(
    ql.SimpleQuote(S)
    )
    flat_ts = ql.YieldTermStructureHandle(
        ql.FlatForward(start, r, day_count)
    )
    dividend_yield = ql.YieldTermStructureHandle(
        ql.FlatForward(start, q, day_count)
    )
    flat_vol_ts = ql.BlackVolTermStructureHandle(
        ql.BlackConstantVol(start, calendar, sigma, day_count)
    )
    bsm_process = ql.BlackScholesMertonProcess(spot_handle, 
                                           dividend_yield, 
                                           flat_ts, 
                                           flat_vol_ts)
    numPaths=10000
    option.setPricingEngine(ql.MCDiscreteArithmeticAPEngine(bsm_process,'lowdiscrepancy',
                                                            antitheticVariate=True,
                                                            controlVariate=True,
                                                            requiredSamples=numPaths,seed=10086))
    return option

def asian_arithmetic_call(S,K,r,q,tau,sigma,freq='1M',**kwargs):
    option=_asian_arithmetic_call(S,K,r,q,tau,sigma,freq)
    epsilon=1e-4
    option_S_upper=_asian_arithmetic_call(S+epsilon,K,r,q,tau,sigma,freq)
    option_S_lower=_asian_arithmetic_call(S-epsilon,K,r,q,tau,sigma,freq)
    option_sigam_upper=_asian_arithmetic_call(S,K,r,q,tau,sigma+epsilon,freq)
    option_sigam_lower=_asian_arithmetic_call(S,K,r,q,tau,sigma-epsilon,freq)
    option_tau_upper=_asian_arithmetic_call(S,K,r,q,tau+100*epsilon,sigma,freq)
    option_tau_lower=_asian_arithmetic_call(S,K,r,q,tau-100*epsilon,sigma,freq)
    option_r_upper=_asian_arithmetic_call(S,K,r+epsilon,q,tau,sigma,freq)
    option_r_lower=_asian_arithmetic_call(S,K,r-epsilon,q,tau,sigma,freq)
    return {
        'PV':option.NPV(),
        'delta':(option_S_upper.NPV()-option_S_lower.NPV())/(2*epsilon),
        'gamma':(option_S_upper.NPV()-2*option.NPV()+option_S_lower.NPV())/(epsilon**2),
        'vega':(option_sigam_upper.NPV()-option_sigam_lower.NPV())/(2*epsilon),
        'theta':(option_tau_upper.NPV()-option_tau_lower.NPV())/(2*100*epsilon),
        'rho':(option_r_upper.NPV()-option_r_lower.NPV())/(2*epsilon),
        }

def _asian_arithmetic_put(S,K,r,q,tau,sigma,freq='1M'):
    # specify timeframe
    tau_day=int(round(tau*365))
    start=ql.Date(1,1,2000)
    ql.Settings.instance().evaluationDate=start
    maturity=start+tau_day
    day_count = ql.Actual365Fixed()
    calendar = ql.UnitedStates()
    period=ql.Period(freq)
    asianFutureFixingDates = [start+period]
    while asianFutureFixingDates[-1]<=maturity:
        asianFutureFixingDates.append(asianFutureFixingDates[-1]+period)
    asianFutureFixingDates.pop(-1)
    # specify option type
    option_type = ql.Option.Put
    payoff = ql.PlainVanillaPayoff(option_type, K)
    exercise = ql.EuropeanExercise(maturity)
    option = ql.DiscreteAveragingAsianOption(ql.Average().Arithmetic,
                                             0, 0, 
                                             asianFutureFixingDates, 
                                             payoff, 
                                             exercise)
    # specify pricing method
    spot_handle = ql.QuoteHandle(
    ql.SimpleQuote(S)
    )
    flat_ts = ql.YieldTermStructureHandle(
        ql.FlatForward(start, r, day_count)
    )
    dividend_yield = ql.YieldTermStructureHandle(
        ql.FlatForward(start, q, day_count)
    )
    flat_vol_ts = ql.BlackVolTermStructureHandle(
        ql.BlackConstantVol(start, calendar, sigma, day_count)
    )
    bsm_process = ql.BlackScholesMertonProcess(spot_handle, 
                                           dividend_yield, 
                                           flat_ts, 
                                           flat_vol_ts)
    numPaths=10000
    option.setPricingEngine(ql.MCDiscreteArithmeticAPEngine(bsm_process,'lowdiscrepancy',
                                                            antitheticVariate=True,
                                                            controlVariate=True,
                                                            requiredSamples=numPaths,seed=10086))
    return option

def asian_arithmetic_put(S,K,r,q,tau,sigma,freq='1M',**kwargs):
    option=_asian_arithmetic_put(S,K,r,q,tau,sigma,freq)
    epsilon=1e-4
    option_S_upper=_asian_arithmetic_put(S+epsilon,K,r,q,tau,sigma,freq)
    option_S_lower=_asian_arithmetic_put(S-epsilon,K,r,q,tau,sigma,freq)
    option_sigam_upper=_asian_arithmetic_put(S,K,r,q,tau,sigma+epsilon,freq)
    option_sigam_lower=_asian_arithmetic_put(S,K,r,q,tau,sigma-epsilon,freq)
    option_tau_upper=_asian_arithmetic_put(S,K,r,q,tau+100*epsilon,sigma,freq)
    option_tau_lower=_asian_arithmetic_put(S,K,r,q,tau-100*epsilon,sigma,freq)
    option_r_upper=_asian_arithmetic_put(S,K,r+epsilon,q,tau,sigma,freq)
    option_r_lower=_asian_arithmetic_put(S,K,r-epsilon,q,tau,sigma,freq)
    return {
        'PV':option.NPV(),
        'delta':(option_S_upper.NPV()-option_S_lower.NPV())/(2*epsilon),
        'gamma':(option_S_upper.NPV()-2*option.NPV()+option_S_lower.NPV())/(epsilon**2),
        'vega':(option_sigam_upper.NPV()-option_sigam_lower.NPV())/(2*epsilon),
        'theta':(option_tau_upper.NPV()-option_tau_lower.NPV())/(2*100*epsilon),
        'rho':(option_r_upper.NPV()-option_r_lower.NPV())/(2*epsilon),
        }

def asian_continuous_geometric_call(S,K,r,q,tau,sigma,**kwargs):
    # specify timeframe
    tau_day=int(round(tau*365))
    start=ql.Date(1,1,2000)
    ql.Settings.instance().evaluationDate=start
    maturity=start+tau_day
    day_count = ql.Actual365Fixed()
    calendar = ql.UnitedStates()
    # specify option type
    option_type = ql.Option.Call
    payoff = ql.PlainVanillaPayoff(option_type, K)
    exercise = ql.EuropeanExercise(maturity)
    option = ql.ContinuousAveragingAsianOption(ql.Average().Geometric,
                                             payoff, 
                                             exercise)
    # specify pricing method
    spot_handle = ql.QuoteHandle(
    ql.SimpleQuote(S)
    )
    flat_ts = ql.YieldTermStructureHandle(
        ql.FlatForward(start, r, day_count)
    )
    dividend_yield = ql.YieldTermStructureHandle(
        ql.FlatForward(start, q, day_count)
    )
    flat_vol_ts = ql.BlackVolTermStructureHandle(
        ql.BlackConstantVol(start, calendar, sigma, day_count)
    )
    bsm_process = ql.BlackScholesMertonProcess(spot_handle, 
                                           dividend_yield, 
                                           flat_ts, 
                                           flat_vol_ts)
    option.setPricingEngine(ql.AnalyticContinuousGeometricAveragePriceAsianEngine(bsm_process))
    return {
    'PV':option.NPV(),
    'delta':option.delta(),
    'gamma':option.gamma(),
    'vega':option.vega(),
    'theta':option.theta(),
    'rho':option.rho(),
    }

def asian_continuous_geometric_put(S,K,r,q,tau,sigma,**kwargs):
    # specify timeframe
    tau_day=int(round(tau*365))
    start=ql.Date(1,1,2000)
    ql.Settings.instance().evaluationDate=start
    maturity=start+tau_day
    day_count = ql.Actual365Fixed()
    calendar = ql.UnitedStates()
    # specify option type
    option_type = ql.Option.Put
    payoff = ql.PlainVanillaPayoff(option_type, K)
    exercise = ql.EuropeanExercise(maturity)
    option = ql.ContinuousAveragingAsianOption(ql.Average().Geometric,
                                             payoff, 
                                             exercise)
    # specify pricing method
    spot_handle = ql.QuoteHandle(
    ql.SimpleQuote(S)
    )
    flat_ts = ql.YieldTermStructureHandle(
        ql.FlatForward(start, r, day_count)
    )
    dividend_yield = ql.YieldTermStructureHandle(
        ql.FlatForward(start, q, day_count)
    )
    flat_vol_ts = ql.BlackVolTermStructureHandle(
        ql.BlackConstantVol(start, calendar, sigma, day_count)
    )
    bsm_process = ql.BlackScholesMertonProcess(spot_handle, 
                                           dividend_yield, 
                                           flat_ts, 
                                           flat_vol_ts)
    option.setPricingEngine(ql.AnalyticContinuousGeometricAveragePriceAsianEngine(bsm_process))
    return {
    'PV':option.NPV(),
    'delta':option.delta(),
    'gamma':option.gamma(),
    'vega':option.vega(),
    'theta':option.theta(),
    'rho':option.rho(),
    }


# ## Visualization

# In[5]:


# =============================================================================
# Visualization
# =============================================================================
def option_plot(option_fcn,changing_param_name,changing_param_bound,other_params,num=100,**kwargs):
    import numpy as np
    import pandas as pd
    import numba as nb
    import plotly.express as px
    from math import ceil
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go
    params=other_params.copy()
    def f(x):
        params[changing_param_name]=x
        return option_fcn(**params)
    f=np.vectorize(f)
    x=np.linspace(changing_param_bound[0],changing_param_bound[1],num)
    y=f(x)
    df=pd.DataFrame(np.array(list(map(lambda x:list(x.values()),y))),columns=y[0].keys(),index=x)
    df.index.name=changing_param_name
#     for i,(name,col) in enumerate(df.iteritems()):
#         fig=px.line(col,y=name,color_discrete_sequence=[px.colors.qualitative.D3[i%10]])
#         fig.show()
    col_num=3
    fig = make_subplots(rows=ceil(len(df.columns)/col_num), cols=col_num)
    for i,(name,col) in enumerate(df.iteritems()):
        fig.add_trace(go.Scatter(x=df.index,y=col),row=i//col_num+1,col=i%col_num+1)
        fig.update_xaxes(title_text=changing_param_name, row=i//col_num+1, col=i%col_num+1)
        fig.update_yaxes(title_text=name, row=i//col_num+1, col=i%col_num+1)
    fig.update_layout(showlegend=False,height=ceil(len(df.columns)/col_num)*300, width=col_num*400)
    fig.show()
    return df


# In[33]:


def interact_plot(changing_param_name,other_params,num=100,**kwargs):
    from functools import partial
    import numpy as np
    from ipywidgets import interact
    from ipywidgets import FloatSlider, Dropdown, FloatRangeSlider
    min_dict={
        'S':0,
        'K':0,
        'X':0,
        'r':-0.05,
        'q':0,
        'tau':0,
        'sigma':0
    }
    max_dict={
        'S':200,
        'K':200,
        'X':200,
        'r':0.2,
        'q':0.5,
        'tau':5,
        'sigma':1
    }
    option_list=[
        ('European Call',vanilla_euro_call),
        ('European Put',vanilla_euro_put),
        ('American Call',vanilla_ame_call),
        ('American Put',vanilla_ame_put),
        ('Barrier UpIn Call',barrier_upin_call),
        ('Barrier UpIn Put',barrier_upin_put),
        ('Barrier UpOut Call',barrier_upout_call),
        ('Barrier UpOut Put',barrier_upout_put),
        ('Barrier DownIn Call',barrier_downin_call),
        ('Barrier DownIn Put',barrier_downin_put),
        ('Barrier DownOut Call',barrier_downout_call),
        ('Barrier DownOut Put',barrier_downout_put),
        ('Asian Geometric Call',asian_geometric_call),
        ('Asian Geometric Put',asian_geometric_put),
        ('Asian Arithmetic Call',asian_arithmetic_call),
        ('Asian Arithmetic Put',asian_arithmetic_put),
        ('Asian Continuous Geometric Call',asian_continuous_geometric_call),
        ('Asian Continuous Geometric Put',asian_continuous_geometric_put)
    ]
    widgets_dict={'option_fcn': Dropdown(
        options=option_list,
        description='Option Type:'
    ),
                  'changing_param_bound':FloatRangeSlider(
                     value=[min_dict[changing_param_name]+(
                         max_dict[changing_param_name]-min_dict[changing_param_name])*0.25,
                           min_dict[changing_param_name]+(
                               max_dict[changing_param_name]-min_dict[changing_param_name])*0.75],
                      min=min_dict[changing_param_name],
                      max=max_dict[changing_param_name],
                      step=10**(round(np.log10(max_dict[changing_param_name]-min_dict[changing_param_name]))-2),
                      description='x range:',
                      disabled=False,
                      continuous_update=False,
                      orientation='horizontal',
                      readout=True,
                      readout_format='.2f' if np.log10(max_dict[changing_param_name]-min_dict[changing_param_name])<1 else 'd',
                 )}
    for key,value in other_params.items():
        widgets_dict[key]=FloatSlider(
                value=value,
                min=min_dict[key],
                max=max_dict[key],
                step=10**(round(np.log10(max_dict[key]-min_dict[key]))-2),
                description=key,
                disabled=False,
                continuous_update=False,
                orientation='horizontal',
                readout=True,
                readout_format='.3f',
        )
    def f(option_fcn,changing_param_bound,**other_params):
        return option_plot(option_fcn,changing_param_name,changing_param_bound,other_params,num=100,**kwargs)
    interact(f,**widgets_dict)


# # Interact Plot

# In[34]:


changing_param_name='S'
changing_param_bound=(50,119)
other_params={
    'S':100,
    'K':100,
    'X':120,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
interact_plot(changing_param_name,other_params)


# # Vanilla European

# ## ...Vs S, call

# In[8]:


changing_param_name='S'
changing_param_bound=(50,150)
other_params={
    'S':100,
    'K':100,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(vanilla_euro_call,changing_param_name,changing_param_bound,other_params)


# ## ...Vs S, put

# In[9]:


changing_param_name='S'
changing_param_bound=(50,150)
other_params={
    'S':100,
    'K':100,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(vanilla_euro_put,changing_param_name,changing_param_bound,other_params)


# ## ...Vs K, call

# In[10]:


changing_param_name='K'
changing_param_bound=(50,150)
other_params={
    'S':100,
    'K':100,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(vanilla_euro_call,changing_param_name,changing_param_bound,other_params)


# ## ...Vs K, put

# In[11]:


changing_param_name='K'
changing_param_bound=(50,150)
other_params={
    'S':100,
    'K':100,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(vanilla_euro_put,changing_param_name,changing_param_bound,other_params)


# ## ...Vs $\sigma$, call

# In[12]:


changing_param_name='sigma'
changing_param_bound=(1e-8,1)
other_params={
    'S':100,
    'K':100,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(vanilla_euro_call,changing_param_name,changing_param_bound,other_params)


# In[13]:


changing_param_name='sigma'
changing_param_bound=(1e-8,1)
other_params={
    'S':150,
    'K':100,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(vanilla_euro_call,changing_param_name,changing_param_bound,other_params)


# ## ...Vs $\sigma$, put

# In[14]:


changing_param_name='sigma'
changing_param_bound=(1e-8,1)
other_params={
    'S':100,
    'K':100,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(vanilla_euro_put,changing_param_name,changing_param_bound,other_params)


# ## ...Vs $\tau$, call

# In[15]:


changing_param_name='tau'
changing_param_bound=(1e-8,5)
other_params={
    'S':100,
    'K':100,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(vanilla_euro_call,changing_param_name,changing_param_bound,other_params)


# ## ...Vs $\tau$, put

# In[16]:


changing_param_name='tau'
changing_param_bound=(1e-8,5)
other_params={
    'S':100,
    'K':100,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(vanilla_euro_put,changing_param_name,changing_param_bound,other_params)


# # Vanilla American

# ## ...Vs S, call

# In[17]:


changing_param_name='S'
changing_param_bound=(50,150)
other_params={
    'S':100,
    'K':100,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(vanilla_ame_call,changing_param_name,changing_param_bound,other_params)


# ## ...Vs S, put

# In[18]:


changing_param_name='S'
changing_param_bound=(50,150)
other_params={
    'S':100,
    'K':100,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(vanilla_ame_put,changing_param_name,changing_param_bound,other_params)


# ## ...Vs K, call

# In[19]:


changing_param_name='K'
changing_param_bound=(50,150)
other_params={
    'S':100,
    'K':100,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(vanilla_ame_call,changing_param_name,changing_param_bound,other_params)


# ## ...Vs K, put

# In[20]:


changing_param_name='K'
changing_param_bound=(50,150)
other_params={
    'S':100,
    'K':100,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(vanilla_ame_put,changing_param_name,changing_param_bound,other_params)


# ## ...Vs $\sigma$, call

# In[21]:


changing_param_name='sigma'
changing_param_bound=(0.001,1)
other_params={
    'S':100,
    'K':100,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(vanilla_ame_call,changing_param_name,changing_param_bound,other_params)


# ## ...Vs $\sigma$, put

# In[22]:


changing_param_name='sigma'
changing_param_bound=(0.001,1)
other_params={
    'S':100,
    'K':100,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(vanilla_ame_put,changing_param_name,changing_param_bound,other_params)


# ## ...Vs $\tau$, call

# In[23]:


changing_param_name='tau'
changing_param_bound=(0.01,1)
other_params={
    'S':100,
    'K':100,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(vanilla_ame_call,changing_param_name,changing_param_bound,other_params)


# ## ...Vs $\tau$, put

# In[24]:


changing_param_name='tau'
changing_param_bound=(0.01,1)
other_params={
    'S':100,
    'K':100,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(vanilla_ame_put,changing_param_name,changing_param_bound,other_params)


# # Barrier

# ## ...Vs S, UpIn call

# In[25]:


changing_param_name='S'
changing_param_bound=(50,119.9)
other_params={
    'S':100,
    'K':100,
    'X':120,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(barrier_upin_call,changing_param_name,changing_param_bound,other_params)


# ## ...Vs S, UpIn put

# In[26]:


changing_param_name='S'
changing_param_bound=(50,119.9)
other_params={
    'S':100,
    'K':100,
    'X':120,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(barrier_upin_put,changing_param_name,changing_param_bound,other_params)


# ## ...Vs S, UpOut call

# In[27]:


changing_param_name='S'
changing_param_bound=(50,119.9)
other_params={
    'S':100,
    'K':100,
    'X':120,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(barrier_upout_call,changing_param_name,changing_param_bound,other_params)


# ## ...Vs S, UpOut put

# In[28]:


changing_param_name='S'
changing_param_bound=(50,119.9)
other_params={
    'S':100,
    'K':100,
    'X':120,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(barrier_upout_put,changing_param_name,changing_param_bound,other_params)


# ## ...Vs S, DownIn call

# In[29]:


changing_param_name='S'
changing_param_bound=(80.1,150)
other_params={
    'S':100,
    'K':100,
    'X':80,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(barrier_downin_call,changing_param_name,changing_param_bound,other_params)


# ## ...Vs S, DownIn put

# In[30]:


changing_param_name='S'
changing_param_bound=(80.1,150)
other_params={
    'S':100,
    'K':100,
    'X':80,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(barrier_downin_put,changing_param_name,changing_param_bound,other_params)


# ## ...Vs S, DownOut call

# In[31]:


changing_param_name='S'
changing_param_bound=(80.1,150)
other_params={
    'S':100,
    'K':100,
    'X':80,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(barrier_downout_call,changing_param_name,changing_param_bound,other_params)


# ## ...Vs S, DownOut put

# In[32]:


changing_param_name='S'
changing_param_bound=(80.1,150)
other_params={
    'S':100,
    'K':100,
    'X':80,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(barrier_downout_put,changing_param_name,changing_param_bound,other_params)


# # Asian

# ## ...Vs S, Geometric call

# In[27]:


changing_param_name='S'
changing_param_bound=(50,150)
other_params={
    'S':100,
    'K':100,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(asian_geometric_call,changing_param_name,changing_param_bound,other_params)


# ## ...Vs S, Geometric put

# In[28]:


changing_param_name='S'
changing_param_bound=(50,150)
other_params={
    'S':100,
    'K':100,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(asian_geometric_put,changing_param_name,changing_param_bound,other_params)


# ## ...Vs S, Arithmetic call

# In[29]:


changing_param_name='S'
changing_param_bound=(50,150)
other_params={
    'S':100,
    'K':100,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(asian_arithmetic_call,changing_param_name,changing_param_bound,other_params)


# ## ...Vs S, Arithmetic put

# In[30]:


changing_param_name='S'
changing_param_bound=(50,150)
other_params={
    'S':100,
    'K':100,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(asian_arithmetic_put,changing_param_name,changing_param_bound,other_params)


# ## ...Vs S, Continuous Geometric call

# In[31]:


changing_param_name='S'
changing_param_bound=(50,150)
other_params={
    'S':100,
    'K':100,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(asian_continuous_geometric_call,changing_param_name,changing_param_bound,other_params)


# ## ...Vs S, Continuous Geometric put

# In[32]:


changing_param_name='S'
changing_param_bound=(50,150)
other_params={
    'S':100,
    'K':100,
    'r':0.01,
    'q':0,
    'tau':1,
    'sigma':0.2
    }
option_plot(asian_continuous_geometric_put,changing_param_name,changing_param_bound,other_params)


# In[ ]:




