module AbstractSDRsSpectrum
# ----------------------------------------------------
## --- Loading Modules
# ----------------------------------------------------
# --- Core Julia modules
using Plots
gr()
using Printf
using FFTW
using AbstractSDRs 
using Statistics
using LinearAlgebra

# --- Package for GUI 
using Blink;        # For GUI in electron 
using Interact;      # For Widget and interactions 
using PlotlyJS; 

# using Suppressor ## To avoid spam when rate is high ?? 
# # export UHD_LOG_FASTPATH_DISABLE=1

# --- Define custom macros
export @updateCarrierFreq
export @updateGain
export @updateBand
export @updateMean
export @updateLim


# ----------------------------------------------------
# --- GUI functions 
# ---------------------------------------------------- 
include("gui.jl");


# ----------------------------------------------------
## --- Processing functions
# ----------------------------------------------------
""" hostSpectrum
---
Acquire signal from USRP and compute the spectrum. Plot the obtained PSD  on a GR plot figure
# --- Syntax
hostSpectrum(nFFT)
# --- Input parameters
- radio         : Radio object
- nFFT			: FFT size
# --- Output parameters
- []
# ---
# v 1.0 - Robin Gerzaguet.
"""
function hostSpectrum(radioAll,nFFT,nbSegMean;maxHold=false);
    # --- Update USRP config.
    # We define global variable to update USRP config
    global doTask	= 1;
    global changed  = 0;
    global newMean  = nbSegMean;
    global newLims  = (-60,20);
    radio = radioAll.rx;
    nbSamples = radio.packetSize;
    global newSamplingRate  = radio.samplingRate;
    # --- Define local parameters
    yLim	        = newLims;
    sF 	            = zeros(Cfloat,nFFT);
    y	            = zeros(Cfloat,nFFT);
    maxH            = zeros(Cfloat,nFFT);
    xAx		        = ((0:nFFT-1) ./ nFFT .-0.5) .*  round(radio.samplingRate,digits=2);
    fMHz            = round(radio.carrierFreq / 1e6, digits=3) ;
    buffer          = zeros(Complex{Cfloat},nbSamples); 
    sig				= zeros(Complex{Cfloat},nFFT); 
    P               = plan_fft(sig;flags=FFTW.PATIENT);
    α 				= 1 / nbSegMean;
    cnt				= 0;
    try 
        while(true)
            # --- Getting signal 
            recv!(buffer,radioAll);
            if nbSamples < nFFT 
                sig[1:nbSamples] .= buffer;
            else 
                sig[1:nFFT] .= buffer[1:nFFT];
            end
            # --- Apply | FFT(.) | ^2
            y  .= abs2.(fftshift(fft(sig)));
            # --- Averaging
            sF .= (1-α) .* sF .+ α .* y; 
            # --- Configuration axis 
            if changed == 1
                nbSegMean   = newMean;
                yLim        = newLims;
                # --- Recalculate x axis based on the new configuration 
                sleep(0.001);
                xAx     = ((0:nFFT-1) ./ nFFT .-0.5) .* round(newSamplingRate,digits=2);
                fMHz	= round(radio.carrierFreq / 1e6, digits=3) ;
                # --- Relock the process until a new change is detected
                changed = 0;
                α 		= 1/nbSegMean;
                sF 		.= 0;
                maxH 	.= 0;
            end
            # --- Update plot 
            plt		= Plots.plot(xAx/1e6,10.0*log10.(sF),title="Spectrum of $(round(radio.samplingRate/1e6,digits=2)) MHz @ $(round(radio.carrierFreq / 1e6, digits=3)) MHz ",xlabel="Frequency [MHz]",ylabel="Power",label="",ylims=yLim,reuse=false);
            if maxHold 
                maxH 	= max.(maxH,sF);
                Plots.plot!(plt,xAx/1e6,10.0*log10.(maxH),label="");
            end
            if mod(cnt,20) == 0
                plt.attr[:size]=(1200,800)
                display(plt);
            end
            cnt+=1;
            yield();
            # --- Interruption manager 
            if doTask != 1
                break;
            end
        end
    catch exception 
        rethrow(exception);
    end
    # --- Release USRP 
    close(radioAll);
end

"""  hostWaterfall
---
Acquire signal from USRP and compute the spectrum. Plot the obtained PSD  on a GR plot figure
# --- Syntax
hostSpectrum(nFFT)
# --- Input parameters
- radio         : Radio object
- nFFT			: FFT size
# --- Output parameters
- []
# ---
# v 1.0 - Robin Gerzaguet.
"""
function hostWaterfall(radio,nFFT,nbSegMean);
    # --- Update USRP config.
    # We define global variable to update USRP config
    global doTask	= 1;
    global changed  = 0;
    global newMean  = nbSegMean;
    global newLims  = (-60,20);
    global newSamplingRate  = radio.samplingRate;
    # --- Define local parameters
    yLim	        = newLims;
    sF 	            = zeros(Cfloat,nFFT,nbSegMean);
    y	            = zeros(Cfloat,nFFT);
    xAx		        = ((0:nFFT-1) ./ nFFT .-0.5) .*  round(radio.samplingRate,digits=2);
    yAy 			= round.((0:nbSegMean-1)./ radio.samplingRate .* 1e3,digits=3)
    fMHz            = round(radio.carrierFreq / 1e6, digits=3) ;
    sig				= zeros(Complex{Cfloat},nFFT); 
    P               = plan_fft(sig;flags=FFTW.PATIENT);
    tSeg 			= 1/ newSamplingRate * nFFT * 1000;
    yAy 	= (0:nbSegMean-1).*tSeg;
    try 
        while(true)
            # --- Getting signal 
            recv!(sig,radio);
            # --- Apply | FFT(.) | ^2
            y  .= abs2.(fftshift(fft(sig)));
            # --- Averaging
            sF	   	= circshift(sF,(0,-1));
            sF[:,end] = y;
            # --- Configuration axis 
            if changed == 1
                # # nbSegMean   = newMean;
                yLim        = newLims;
                # --- Recalculate x axis based on the new configuration 
                xAx     = ((0:nFFT-1) ./ nFFT .-0.5) .* round(newSamplingRate,digits=2);
                # yAy 	= (0:nbSegMean*nFFT-1) ./ radio.samplingRate .* 1e3;
                fMHz	= round(radio.carrierFreq / 1e6, digits=3) ;
                tSeg 	= 1/newSamplingRate * nFFT*1000;
                yAy 	= (0:nbSegMean-1).*tSeg;
                # --- Relock the process until a new change is detected
                changed = 0;
                sF 		.= 0;
            end
            # --- Update plot 
            # plt		= heatmap(xAx/1e6,(0:nbSegMean-1),10.0*log10.(sF)',title="Spectrum of $(round(radio.samplingRate/1e6,digits=2)) MHz @ $(round(radio.carrierFreq / 1e6, digits=3)) MHz ",xlabel="Frequency [MHz]",ylabel="Time [burst]",clim=yLim,label="",reuse=false);
            plt		= heatmap(xAx/1e6,yAy,10.0*log10.(sF)',title="Spectrum of $(round(radio.samplingRate/1e6,digits=2)) MHz @ $(round(radio.carrierFreq / 1e6, digits=3)) MHz ",xlabel="Frequency [MHz]",ylabel="Time [ms]",clim=yLim,label="",reuse=false);

            # plt = heatmap(transpose(sF));
            # end
            plt.attr[:size]=(1200,800)
            gui(plt);
            yield();
            # gui(plt);
            # --- Interruption manager 
            if doTask != 1
                break;
            end
        end
    catch exception ;
        rethrow(exception);
    end
    # --- Release USRP 
    close(radio);
end

function processing!(out,sig,internal,P)
    # --- Plan FFT 
    mul!(internal,P,sig);
    # Shift 
    sig .= fftshift(internal)
    # --- |.|^2 
    out .= abs2.(sig);
end



""" killSpectrum
---
Interrupt hostSpectrum function when run in asynchronous mode
# --- Syntax
killSpectrum()
# --- Input parameters
- []
# --- Output parameters
- []
# ---
# v 1.0 - Robin Gerzaguet.
"""
function stop()
    global doTask = 0;
end
function killSpectrum()
    global doTask			  = 0;
end

""" @updateCarrierFreq
---
Update carrier frequency during T/F analysis
# --- Syntax
@updateCarrierFreq xxx
# --- Input parameters
- xxx	: Value of new desired carrier frequency
# --- Output parameters
- []
# ---
# v 1.0 - Robin Gerzaguet.
"""
macro updateCarrierFreq(param)
    # --- Updating carrier frequency
    global changed = 1;
    global newCarrierFreq = param;
    # --- Calling routine to update radio
    updateCarrierFreq!(radio,param);
end

""" @updateGain
---
Macro to dynamically update Rx gain during T/F grid
# --- Syntax
@updateGain xxx
# --- Input parameters
- xxx	: New value of Rx gain
# --- Output parameters
- []
# ---
# v 1.0 - Robin Gerzaguet.
"""
macro updateGain(param)
    # --- Updating carrier frequency
    global newGain = param;
    global changed = 1;
    # --- Calling routine to update radio
    updateGain!(radio,param);
end


""" @updateBand
---
Dynamically update Rx sample rate
# --- Syntax
@updateBand xxx
# --- Input parameters
- xxx	: New desired value for sample rate
# --- Output parameters
- []
# ---
# v 1.0 - Robin Gerzaguet.
"""
macro updateBand(param)
    # --- Updating carrier frequency
    global newSamplingRate = param;
    global changed = 1;
    # --- Calling routine to update radio
    updateSamplingRate!(radio,newSamplingRate);
end

""" @updateLim
---
Update the y axis limit view dynamically
# --- Syntax
@updateLim (xxx1,xxx2)
# --- Input parameters
- xxx1	  : New minimal y value
- xxx2	  : New maximal  y value
# --- Output parameters
- []
# ---
# v 1.0 - Robin Gerzaguet.
"""
macro updateLim(param)
    global changed = 1;
    global newLims=eval(param);
end

"""	@updateMean
---
Update the averaging system to compute spectrum
# --- Syntax
@updateMean xxx
# --- Input parameters
- xxx	: New number of mean
# --- Output parameters
- []
# ---
# v 1.0 - Robin Gerzaguet.
"""
macro updateMean(param)
    global changed = 1;
    global newMean = param;
end


# ----------------------------------------------------
## --- Main call
# ----------------------------------------------------
""" main
---
Main routine to display Power spectral density
# --- Syntax
main()
# --- Input parameters
- []
# --- Output parameters
- []
# ---
# v 1.0 - Robin Gerzaguet.
"""
function start(sdr;kwargs...)
    main(sdr;kwargs...);
end
function main(sdr::Union{Symbol,String};carrierFreq=770e6,gain=20,samplingRate=80e6,maxHold::Bool=false,waterfall::Bool=false,kwargs...)
    @info "Spectrum display -- "
    # --- Simulation parameters
    global newCarrierFreq	  = carrierFreq;		  # --- Starting frequency
    global newGain		      = gain;			  # --- Analog Rx gain
    global newSamplingRate	  = samplingRate;		  # --- Sampling
    nFFT 		              = 1024;
    # --- Update radio configuration
    global radio			= openSDR(sdr,newCarrierFreq,newSamplingRate,newGain;kwargs...); 
    if !waterfall
        @async hostSpectrum(radio,nFFT,32;maxHold=maxHold);
    else 
        @async hostWaterfall(radio,nFFT,32);
    end
end





end # module
