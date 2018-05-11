%Simulation of an epidemic model
%inputs : mod (model structure), param (parameter structure)
%output : modsol (model solution structure)
function modsol = simulatemodel(mod,param)

    nc = numel(mod.comp); %number of model compartments
    ns = numel(mod.sink); %number of model sinks
    na = numel(mod.accum); %number of model accumulators
    xi = zeros(nc + ns + na,1); %initial state vector
    
    contactsplit = containers.Map;
    contactkeys = keys(mod.contact);
        
    for mm = 1 : numel(contactkeys)
        testkeym = char(contactkeys(mm));
        vecm = strsplit(testkeym,{'-',','});
        jm = find(strcmp(mod.comp,vecm(1)));
        for sm = 2:length(vecm)
            contactsplit(char(strcat(vecm(1),'-',vecm(sm)))) = mod.contact(testkeym);
        end
    end
    
    compkeys = keys(mod.init.comp);
    for k = 1:numel(compkeys)
        testkey = char(compkeys(k));
        c = find(strcmp(mod.comp,testkey));
        xi(c) = mod.init.comp(testkey);
    end
    
    sinkkeys = keys(mod.init.sink);
    for k = 1:numel(sinkkeys);
        testkey = char(sinkkeys(k));
        c = find(strcmp(mod.sink,testkey));
        xi(nc + c) = mod.init.sink(testkey);
    end
    
    accumkeys = keys(mod.init.accum);
    for k = 1:numel(accumkeys);
        testkey = char(accumkeys(k));
        c = find(strcmp(mod.accum,testkey));
        xi(nc + ns  + c) = mod.init.accum(testkey);
    end
    
    %determine number of simulation periods
    ntime = 1;
    for k = 1 : numel(mod.event)-1
        tspan = [mod.event(k).time mod.event(k+1).time];
        tsol = tspan(1):mod.timestep:tspan(2);
        ntime = ntime + numel(tsol)-1;
    end
    
    fullsol = zeros(nc+ns+na, ntime);
    fullsol(:,1) = xi;
        
    xo = xi;
    timecounter = 1;
    tmod = mod;
    tevents = mod.event;

    %first simulation event
    xo = callEvent(tmod,param,tevents,1,xo,nc,ns,na);
    
    %main simulation loop
    for k = 1 : numel(tevents)-1
        tspan = [tevents(k).time tevents(k+1).time]
        tsol = tspan(1):tmod.timestep:tspan(2);
        sol = ode45(@modeldiff,tspan,xo);
        xsol = deval(sol,tsol);       
        fullsol(:,timecounter:timecounter+numel(tsol)-1) = xsol(:,1:end);
        timecounter = timecounter + numel(tsol)-1;   
        xo = xsol(:,end);       
        [xo tmodout] = callEvent(tmod,param,tevents,k+1,xo,nc,ns,na);
    end
    
    %empty placeholder for solution
    modsol.comp = containers.Map;
    modsol.sink = containers.Map;
    modsol.accum = containers.Map;
    
    for k = 1:numel(tmod.comp)
        modsol.comp(char(tmod.comp(k))) = fullsol(k,:);
    end
    
    for k = 1:numel(tmod.sink)
        modsol.sink(char(tmod.sink(k))) = fullsol(nc+k,:);
    end
    
    for k = 1:numel(tmod.accum)
        modsol.accum(char(tmod.accum(k))) = fullsol(nc+ns+k,:);
    end
    
    modsol.time = tmod.event(1).time:tmod.timestep:tmod.event(numel(tmod.event)).time;
    
    %subfunction to compute derivative of state vector
    function dx = modeldiff(t,x)

        dx = zeros(nc+ns+na,1);
        
        %sources
        inkeys = keys(tmod.enter);
        
        for m = 1 : numel(inkeys)
            testkey = char(inkeys(m));
            vec = strsplit(testkey,'-');
            if(length(vec) == 1)
                %constant recruitment
                s = find(strcmp(tmod.comp, vec(1)));
                dx(s) = dx(s) + tmod.enter(testkey);
            elseif(length(vec) == 2)
                %linear recruitment
                j = find(strcmp(tmod.comp, vec(1)));
                s = find(strcmp(tmod.comp, vec(2)));
                dx(s) = dx(s) + tmod.enter(testkey)*x(j);
            else
                error('Error in source keys.');
            end
        end
        
        
        %sinks
        outkeys = keys(tmod.exit);
        
        for m = 1 : numel(outkeys)
            testkey = char(outkeys(m));
            vec = strsplit(testkey,'-');
            j = find(strcmp(tmod.comp,vec(1)));
            s = find(strcmp(tmod.sink,vec(2)));
            temp = tmod.exit(testkey) * x(j);
            dx(j) = dx(j) - temp;
            dx(nc + s) = dx(nc + s) + temp;
        end
        
        
        %links
        linkkeys = keys(tmod.flow);
        
        for m = 1 : numel(linkkeys)
            testkey = char(linkkeys(m));
            vec = strsplit(testkey,'-');
            j = find(strcmp(tmod.comp,vec(1)));
            s = find(strcmp(tmod.comp,vec(2)));
            testlinkval = tmod.flow(testkey);
            
            if(ischar(testlinkval))
                fh = str2fun(testlinkval);
                linkval = fh(t);
            elseif(isfloat(testlinkval))
                linkval = testlinkval;
            else
                error('Error in link.') 
            end
            
            temp = linkval*x(j);
            dx(j) = dx(j) - temp;
            dx(s) = dx(s) + temp;
            
            %accumulate transitions
            temp2 = find(strcmp(tmod.accum,testkey),1);
            if(~isempty(temp2))
                a = temp2;
                dx(nc + ns + a) = dx(nc + ns + a) + temp;
            end
        end
        
        %new infections
        sizepool = zeros(nc,1);
        contactsplitkeys = keys(contactsplit);
        
        for m = 1 : numel(contactsplitkeys)
            testkey = char(contactsplitkeys(m));
            vec = strsplit(testkey,'-');
            j = find(strcmp(tmod.comp,vec(1)));
            s = find(strcmp(tmod.comp,vec(2)));
            sizepool(j) = sizepool(j) + x(s);
        end
        
        forcekeys = keys(tmod.force);
        
        for m = 1 : numel(forcekeys)
            testkey = char(forcekeys(m));
            vec = strsplit(testkey,{',','-'});
            a = find(strcmp(tmod.comp,vec(1)));
            j = find(strcmp(tmod.comp,vec(2)));
            s = find(strcmp(tmod.comp,vec(3)));
            
            %protection vs incoming infection
            ieff = 0.0;
            testkey2 = char(vec(2));
            if(isKey(tmod.inceff,testkey2))
                testinceff = tmod.inceff(testkey2);
                
                if(ischar(testinceff))
                    fh = str2fun(testinceff);
                    ieff = fh(t);
                elseif(isfloat(testinceff))
                    ieff = testinceff;
                else
                    error('Error in incoming infection efficacy.') 
                end
            end
            
            %prevention of outgoing infection
            oeff = 0.0;
            testkey3 = char(vec(1));
            
            if(isKey(tmod.outeff,testkey3))
                testouteff = tmod.outeff(testkey3);
                
                if(ischar(testouteff))
                    fh = str2fun(testouteff);
                    oeff = fh(t);
                elseif(isfloat(testouteff))
                    oeff = testouteff;
                else
                    error('Error in outgoing infection efficacy.') 
                end
            end
            
            ctact = 0.0;
            testkey4 = char(strcat(vec(2),'-',vec(1)));
            
            if(isKey(contactsplit,testkey4))
                ctact = contactsplit(testkey4);
            end
            
            temp = ctact*(1-ieff)*(1-oeff)*tmod.force(testkey)*x(j)*x(a)/sizepool(j);
            dx(j) = dx(j) - temp;
            dx(s) = dx(s) + temp;
            
            %accumulate infections
            testkey5 = strcat(vec(2),'-',vec(3));
            temp2 = find(strcmp(tmod.accum,testkey5),1);
            
            if(~isempty(temp2))
                a = temp2;
                dx(nc + ns + a) = dx(nc + ns + a) + temp;
            end        
        end   
    end

    %subfunction to call simulation event
    function [xout tmodout] = callEvent( tmodin,param,events, nev,xin,nc,ns,na)

        tmodt = tmodin;
        tmodt.current.name = events(nev).name;
        tmodt.current.time = events(nev).time;
        tmodt.current.comp = containers.Map;
        tmodt.current.sink = containers.Map;
        tmodt.current.accum = containers.Map;
        
        for n = 1 : nc
            tmodt.current.comp(char(tmodt.comp(n))) = xin(n);
        end
        
        for n = 1 : ns
            tmodt.current.sink(char(tmodt.sink(n))) = xin(nc+n);
        end
        
        for n = 1 : na
            tmodt.current.accum(char(tmodt.accum(n))) = xin(nc+ns+n);
        end

        tfc = events(nev).fn;

        if(ischar(tfc))
            fs = str2func(tfc);
            tmodt = fs(tmodt,param);
        elseif(iscell(tfc))

            for nev = 1 : numel(tfc)
                fs = str2func(char(tfc(nev)));
                tmodt = fs(tmodt,param);
            end
        else
            error('Invalid event function.');
        end
        
        xout = xin;
        
        for n = 1 : nc
            xout(n) = tmodt.current.comp(char(tmodt.comp(n)));
        end
        
        for n = 1 : ns
            xout(nc+n) = tmodt.current.sink(char(tmodt.sink(n)));
        end
        
        for n = 1 : na
            xout(nc+ns+n) = tmodt.current.accum(char(tmodt.accum(n)));
        end
        
        tmodt = rmfield(tmodt,'current');
        tmodt.event = events;
        tmodout = tmodt;
    end
end

