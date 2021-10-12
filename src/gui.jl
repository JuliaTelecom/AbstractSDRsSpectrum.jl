# ----------------------------------------------------
# --- Global functions 
# ---------------------------------------------------- 
FLAG    = false
UPDATE  = false
CURRENT = true
MAX_HOLD = false
nbSegMean = 24
radio = nothing

WINDOW_SIZE = (800,700); # Change it by @eval AbstractSDRsSpectrum WINDOW_SIZE=(400,200)
PLOT_SIZE = (700,500);
MAX_BANDWIDTH = 16;

# ----------------------------------------------------
# --- Main GUI core processing function
# ---------------------------------------------------- 
function mainGUI(sdr,carrierFreq,samplingRate,gain=20,nbSegMeanInit=24,p=nothing,maxHold=false;kwargs...)
    # --- Simulation parameters
    global FLAG;
    global UPDATE;
    global nbSegMean = nbSegMeanInit;
    global newCarrierFreq	  = carrierFreq;		# --- Starting frequency
    # --- Duration and buffer size 
    duration                  = 0.5;
    nbSamples                 = Int( duration * samplingRate);
    # --- Update radio configuration
    global radio			= openSDR(sdr,newCarrierFreq,samplingRate,gain;kwargs...); 
    # --- Create FIR filter 
    nFFT      = 1024;
    nbSamples = 1024;
    buffer    = zeros(Complex{Cfloat},nbSamples);
    sig       = zeros(Complex{Cfloat},nFFT);
    sfM         = zeros(Cfloat,nFFT);
    maxS        = -120*ones(Cfloat,nFFT);
    y= zeros(Cfloat,nFFT);
    sF        = zeros(Cfloat,nFFT);
    α         = 1/ nbSegMean;
    x         = ((0:nFFT-1)/nFFT .- 0.5) .* samplingRate ./1e6;
    cnt = 0;
    FLAG = true;
    # --- Audio rendering 
    try 
        while(true) 
            if UPDATE 
                # --- We have update something, need to redo sme stuff 
                nbSegMean == 0 ? α = 1 :  α         = 1/ nbSegMean;
                sfM .= 0;
                x         = ((0:nFFT-1)/nFFT .- 0.5) .* getSamplingRate(radio) ./1e6;
                UPDATE = false;
                maxS .= -120;
                maxHold = MAX_HOLD;
                deletetraces!(p,1);
            end
            # --- We get buffer 
            # --- Getting signal 
           cnt +=1 
            recv!(buffer,radio);
            if nbSamples < nFFT 
                sig[1:nbSamples] .= buffer;
            else 
                sig[1:nFFT] .= buffer[1:nFFT];
            end
            global tmp = sig;
            # --- Apply | FFT(.) | ^2
            sF  .= abs2.(fftshift(fft(sig)));
            # --- Averaging
            sfM .= (1-α) .* sfM .+ α .* sF; 
            y .= 10*log10.(sfM);
            maxS .= max.(maxS,y);
            # --- Plot 
            if mod(cnt,100) == 0
                deletetraces!(p,0);
                if maxHold == true  && cnt > 1
                    deletetraces!(p,1);
                end 
                pl1	  = PlotlyJS.scatter(; x,y, name=" ");
                addtraces!(p,pl1);
                if maxHold == true 
                    maxS .= max.(maxS,y);
                    pl2	  = PlotlyJS.scatter(; x,y=maxS, name=" ");
                    addtraces!(p,pl2);
                    sleep(0.1);
                else 
                    sleep(0.05);
                end 
            end
            yield();
            # Exit 
            if FLAG == false;
                break;
            end
        end
    catch exception;
        close(radio);
        rethrow(exception);
    end 
    close(radio);
end

# ----------------------------------------------------
# --- GUI
# ---------------------------------------------------- 
function gui(;kw1...)
    global FLAG = false;
    global UPDATE = false;
    # ----------------------------------------------------
    # --- Define widgets
    # ---------------------------------------------------- 
    # --- Creating widget for FM frequencies 
    # frequencies = widget(100:4:5000, label="FM frequencies (MHz)");
    widgetFrequencies = textbox(hint=""; value="800",label="Frequency (MHz)  ")
    # widgetFrequencies= widget(100:10:5000, label="Radio Frequencies (MHz)");
    widgetGain        = widget(0:50, label="Radio Gain");
    widgetAverage     = widget(0:512, label="Averaging window");
    widgetBandwidth   = widget(0:1:MAX_BANDWIDTH, label="Bandwidth (MHz)");
    # --- Creating some stuff related to args 
    widgetArgs = textbox(hint=""; value="arg=\"addr=192.168.10.13\"");
    widgetMaxHold = checkbox(false,label="Max Hold");
    # --- Widget for radio configuration 
    sdrList = AbstractSDRs.getSupportedSDRs();
    # options = Observable([:e310,:uhd,:radiosim])
    options = Observable(sdrList);
    wdg = dropdown(options);
    # --- Start button 
    startButton = button("Start !"; value=0)
    configButton= button("Radio Config"; value=0,style = Dict("backgroundColor" => "#636EFA") );
    stopButton = button(" ! Stop"; value=0,style = Dict("backgroundColor" => "#EF553B") );
    # --- 
    (wW,hW) = WINDOW_SIZE;
    # w = Window(Blink.@d(:width=>wW, :height=>hW));
    w = Window();
    title(w,"Spectrum viewer in Julia");
    # --- Horizontal axis 
    layout = hbox(pad(1em,wdg),pad(1em,widgetArgs));
    layout = vbox(layout,widgetGain);
    layout = vbox(layout,widgetAverage);
    layout = vbox(layout,widgetFrequencies);
    layout = vbox(layout,widgetBandwidth);
    layout = hbox(layout,pad(1em,vbox(startButton,configButton,stopButton,widgetMaxHold)));
    # 
    (wW,hW) = PLOT_SIZE;
    plotLayout = Layout(;
                        width = wW,
                        height = hW,
                        title=" Spectrum ",
                        xaxis_title=" Frequency MHz ",
                        yaxis_title=" Power ",
                        xaxis_showgrid=true, yaxis_showgrid=true,
                        yaxis_range =[-80,0]
                        )
    N = 100;
    x = 0:N-1;
    y = zeros(N);
    pl1	  = PlotlyJS.scatter(; x,y, name=" ");
    p = PlotlyJS.plot(pl1,plotLayout)
    layout = vbox(layout,p);
    body!(w,layout);
    # ----------------------------------------------------
    # --- 
    # ---------------------------------------------------- 
    # --- Update Frequency 
    on(updateC,widgetFrequencies);
    # --- Update Gain 
    on(updateG,widgetGain);
    # --- Update average 
    on(updateM,widgetAverage);
    on(updateB,widgetBandwidth);
    v = on(configButton) do click
        c = Dates.format(now(), "HH:MM:SS")  
        println("--------\n $c New config radio");
        print(radio);
    end
    m = on(widgetMaxHold) do mm 
        global MAX_HOLD = mm;
        UPDATE = true;
    end
    s = on(stopButton) do click 
        if FLAG == true 
            startButton["is-loading"][]=false;
            global FLAG = false;
            # close(w);
        end
    end
    # --- Starting routine 
    h = on(startButton) do click 
        startButton["is-loading"][]=true
        #@show click
        #if click == 1 
        #    map!(fromStartToStop,startButton,12,16)
        #end
        sdr =  wdg[];
        carrierFreq = parse(Float64,widgetFrequencies[])*1e6
        # carrierFreq = widgetFrequencies[]*1e6;
        gain        = widgetGain[];
        nbSegMean   = widgetAverage[];
        samplingRate = widgetBandwidth[] *1e6;
        maxHold = widgetMaxHold[];
        # kwargs = (;args=widgetArgs[]);
        # --- Creating a dict with all keywords, given as a string
        str = widgetArgs[] 
        argS = split(str,",") 
        kwargs = Dict()
        for s in argS 
            (k,p_) = split(s,"=")
            p_ = replace(p_,"\""=>"")
            push!(kwargs,Symbol(k)=>p_)
        end
        # kwargs = (;args=widgetArgs[]);
        @async mainGUI(sdr,carrierFreq,samplingRate,gain,nbSegMean,p,maxHold;kwargs...,kw1...);
    end
    success(w.shell.proc)
    global FLAG=false; 
    close(radio);
    # We get the hand back and can close the radio 
    return nothing
end

# ----------------------------------------------------
# --- Callbacks
# ---------------------------------------------------- 
function updateC(val)
    global FLAG, UPDATE, CURRENT
    if FLAG == true 
        if val != ""
            carrierFreq = parse(Float64,val)*1e6;
            if CURRENT == true 
                # When we write new value it can lead to several calls to updateCarrierFreq, not necessary well handle on the SDR side... We implement a soft lock here to be sure we do not call sevearl times updateCarrierFreq
                updateCarrierFreq!(radio,carrierFreq);
                UPDATE = true
                CURRENT = false 
            end
            CURRENT = true
        end
    end
end
function updateG(val)
    global FLAG, UPDATE;
    if FLAG == true 
        gain = val;
        updateGain!(radio,gain);
        UPDATE = true;
    end
end

function updateM(val)
    global nbSegMean = val;
    global UPDATE = true;
end
function updateB(val)
    global FLAG, UPDATE;
    if FLAG == true 
        samplingRate = val*1e6;
        # carrierFreq = val*1e6;
        updateSamplingRate!(radio,samplingRate);
        UPDATE = true;
    end
end



