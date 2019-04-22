import math
import pint
uReg = pint.UnitRegistry()


inf = float("inf")
Wm2 = uReg.W / uReg.m**2
Jm2 = uReg.J / uReg.m**2
mps = uReg.m/uReg.sec

#typical pulse duration in seconds
laserTypes={
	"D": {"desc":"Continuous wave", "typical pulse duration":(0.2*uReg.s, inf), "count":1, },
	"I": {"desc":"Pulsed", "typical pulse duration": (1e-6*uReg.s, 0.2*uReg.s), "count":100, },
	"R": {"desc":"Giant-pulsed", "typical pulse duration": (1e-9*uReg.s, 1e-6*uReg.s), "count":100, },
	"M": {"desc":"Mode-coupled", "typical pulse duration": (0*uReg.s, 1e-9*uReg.s), "count":100, },
}

EN166Robustness={
	"S": "increased",
	"A": 190*mps,
	"B": 120*mps,
	"F": 45*mps,
	"S": "Temperature",
}

manufacturerIDs={
	"RZ": "Росомз",
	"YL": "Yamamoto Kogaku Co.,Ltd",
}


visibleRange = (400*uReg.nm, 1400*uReg.nm)

#nm
EN207={
	(180*uReg.nm, 315*uReg.nm): {"D": 1e-3*Wm2, "IR": 3e1*Jm2, "M":3e10*Jm2, },
	(315*uReg.nm, 1400*uReg.nm): {"D": 1e1*Wm2, "IR": 5e-3*Jm2, "M":1.5e-4*Jm2, },
	(1400*uReg.nm, 1000000*uReg.nm): {"D": 1e3*Wm2, "IR": 1e-4*Jm2, "M":1e11*Jm2, }
}

EN208={ # EN208 is in Watts / Joules!
	False: 1e-3*uReg.W, # CW or pulsed with long pulses
	True: 2e-7*uReg.J # pulsed with short pulses
}
EN208ThresholdPulseDurationContinious = 2e-4 * uReg.s
EN208ThresholdPulseDurationForbidden = 1e-9 * uReg.s

ranges=EN207

pupilDiameter = 7 * uReg.mm
blikReflexTime = 0.25 * uReg.s

def EN60825MPE(t):
	return _EN60825MPE(t.to(uReg.s).magnitude()) * Jm2

def _EN60825MPE(t):
	return 18 * (t**0.75) #Jm-2

def NOHD_(P, D, divergence):
	return math.sqrt(4*P/(math.pi*MPE)-D)/divergence

def computeTransmittance(L):
	return 10**(-L)

def computeAttenuation(L):
	return 10**L

def getRange(wavelength):
	for r, rangeContents in ranges.items():
		if r[0] <= wavelength < r[1]:
			return r, rangeContents


def irradianceToLevel(irradiance, allowedIrradiance):
	return int(math.ceil(math.log10(irradiance/allowedIrradiance)))

def level2Irradiace(level, allowedIrradiance):
	return allowedIrradiance*computeAttenuation(level)

timeInComputation = 5 * uReg.s

materialCorrections={
	"glass": 1.1693,
	"plastic": 1.2233
}

def materialCorrection_(d, material):
	if d > 1.:
		if material in materialCorrections:
			return d**materialCorrections[material]
	
	return 1.

def materialCorrection(d, material):
	return materialCorrection_(d.to(uReg.mm).magnitude, material)

def computeGoogles__(powerOrItsDensity, impulseEnergyOrItsDensity, laserType, impulsePower, isVisible, rangeContents, CWType="D"):
	level = irradianceToLevel(powerOrItsDensity, rangeContents[CWType])
	if impulseEnergyOrItsDensity is not None:
		correction = 1.
		if isVisible:
			correction = (pulseRate*timeInComputation)**0.25
		impulseEnergyOrItsDensity*=correction
		levelPulsed = irradianceToLevel(impulseEnergyOrItsDensity, rangeContents[laserType])
	
		level = max(level, levelPulsed)
	return level

def computeGoogles_(wavelength, laserType, waist, averagePower, impulseEnergy=None, pulseTime=0, pulseRate=0., needView=True):
	r=waist/2
	area=math.pi*r**2
	powerDensity = averagePower / area

	isVisible = visibleRange[0] <= wavelength < visibleRange[1]
	isContinious = laserType == "D"

	EN208Applies = isVisible and needView and (isContinious or (pulseTime > EN208ThresholdPulseDurationForbidden))

	res=[]
	if EN208Applies:
		laserTypeEN208 = (isContinious or pulseTime > EN208ThresholdPulseDurationContinious)
		level = computeGoogles__(averagePower, impulseEnergy, laserTypeEN208, isContinious, isVisible, EN208, CWType=False)
		res.append(encodeFilterEN208(wavelength, level))

	spectralRange, rangeContents = getRange(wavelength)
	level = computeGoogles__(powerDensity, (impulseEnergy / area) if (impulseEnergy is not None) else None, laserType, isContinious, isVisible, rangeContents, CWType="D")
	
	res.append(encodeFilterEN207(wavelength, laserType, level))
	return res


def computeGooglesPulsed(wavelength, laserType, waist, impulseEnergy, pulseTime, pulseRate, needView=True):
	return computeGoogles_(wavelength, laserType, waist, impulseEnergy*pulseRate, impulseEnergy, pulseTime, pulseRate, needView)


def computeGooglesContinious(wavelength, waist, averagePower, needView=True):
	return computeGoogles_(wavelength, "D", waist, averagePower, None, None, None, needView)


def encodeFilterRest(prefix, manufCode, DINCompliance=False, conformiteEuropeene=False, robustness=False):
	return " ".join(filter(None, (prefix, manufCode, "DIN" if DINCompliance else None, "CE" if conformiteEuropeene else None, robustness if robustness else None)))

def encodeFilterEN207_(wavelength:int, laserType:str, L:int, manufCode:str="", DINCompliance=False, conformiteEuropeene=False, robustness=False):
	return encodeFilterRest(" ".join((str(wavelength), laserType, "LB" + str(L))), manufCode, DINCompliance, conformiteEuropeene, robustness)


def encodeFilterEN207(wavelength:int, laserType:str, L:int, manufCode:str="", DINCompliance=False, conformiteEuropeene=False, robustness=False):
	return encodeFilterEN207_(int(round(wavelength.to(uReg.nm).magnitude)), laserType, L, manufCode, DINCompliance, conformiteEuropeene, robustness)

def encodeFilterEN208_(maxPower, maxEnergy, wavelength:int, RB:int, manufCode:str="", DINCompliance=False, conformiteEuropeene=False, robustness=False):
	return encodeFilterRest(" ".join((str(maxPower)+"W", str(maxEnergy)+"J", str(wavelength), "RB"+str(RB))), manufCode, DINCompliance, conformiteEuropeene, robustness)

def encodeFilterEN208(wavelength:int, RB:int, manufCode:str="", DINCompliance=False, conformiteEuropeene=False, robustness=None):
	maxPower=int(math.ceil(level2Irradiace(RB, EN208[False]).to(uReg.W).magnitude))
	maxEnergy=int(math.ceil(level2Irradiace(RB, EN208[True]).to(uReg.J).magnitude))
	return encodeFilterEN208_(
		maxPower,
		maxEnergy,
		int(round(wavelength.to(uReg.nm).magnitude)),
		RB,
		manufCode, DINCompliance, conformiteEuropeene, robustness,
	)
