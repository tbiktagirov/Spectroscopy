function H = integration()

	file = 'E6P1_532'
	delta = 2.9e-1;
	gamma = 1.7e-4;	
	
	%file = 'E6P1_633'
	%delta = 2.9e-1;
	%gamma = 2.1e-4;	
	
	%file = 'P_532'
	%delta = 3.1e-1;
	%gamma = 1.65e-4;

	%file = 'P_633'
	%delta = 3.1e-1;
	%gamma = 2.0e-4;	
	
	%file = 'E6P2_532'
	%delta = 3.5e-1;
	%gamma = 2.4e-4;
	
	%file = 'E6P2_633'
	%delta = 3.7e-1;
	%gamma = 2.9e-4;		

	%file = 'E6_as_grown'
	%delta = 0.18;
	%gamma = 1e-6;	

	
	
	%transform experimental data
	data=load([file '.csv']);
	data(:,1) = data(:,1) + 1041.6893;
	y=data(:,2);
	wavelength=data(:,1);
	lambda0 = 1041.6893;
	w0 = 1e7/lambda0;
	freq=[];
	for i=1:size(wavelength,1)
		freq(i) = 9599.9 + (1/wavelength(i) - 1/lambda0)*1e7; %cm-1
	end
	freq1 = freq * 1.23981e-4 * 1e3; %meV
	
	
	
	%interaction parameters
	B = -1.23;
	C = -0.69;
	c11 = 1076;
	c12 = 125;
	c44 = 576;
	%the following is from Davies 1979 J Phys C :
	v1 = (c11-c12)*B + c44*C; %in meV
	v2 = sqrt(2)*(c11-c12)*B - c44*C/sqrt(2); %in meV


	
	x0 = 1190.15;
	freq=linspace(-4,4,61)+x0;
	H=zeros(size(freq,1),1);
	
	%Monte Carlo integration
	N = 1e5
	for i=1:size(freq,2)
		tmp = 0;
		x=freq(i);
		for j = 1:N
			t1 = rand(1,1);
			t2 = rand(1,1);
			t3 = rand(1,1);
			t4 = rand(1,1);
			tmp = tmp + (lineshape(t1, t2, t3, t4, x, x0, v1, v2, delta, gamma));
		end
		H(i)=tmp/N;
	end
	plot(freq,H/max(H), freq1, y/max(y))
	
	csvwrite(['sim_' file '.csv'], [freq',H'])
	csvwrite(['exp_' file '.csv'], [freq1',y])
	%R = H/max(H) - y/max(y);
