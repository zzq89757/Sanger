function convertBasePosTraceToPerBpTrace(chromData) {
	const { basePos, aTrace } = chromData;
	const traceLength = aTrace.length;
	let startPos = 0;
	let nextBasePos = basePos[1];
	let endPos;
	function setEndPos() {
		if (nextBasePos) {
			endPos = startPos + Math.ceil((nextBasePos - startPos) / 2);
		} else {
			endPos = traceLength;
		}
	}
	setEndPos();
	const baseTraces = [];
	for (let i = 0; i < basePos.length; i++) {
		const tracesForType = {
			aTrace: [],
			tTrace: [],
			gTrace: [],
			cTrace: [],
		};
		baseTraces[i] = tracesForType;
		[
			"aTrace",
			"tTrace",
			"gTrace",
			"cTrace",
			// eslint-disable-next-line no-loop-func
		].forEach((type) => {
			const traceForType = tracesForType[type];
			const traceData = chromData[type];
			for (let j = startPos; j < endPos + correctionAmount; j++) {
				traceForType.push(traceData[j] || 0);
			}
		});
		if (i !== basePos.length - 1) {
			startPos = endPos + correctionAmount;
			nextBasePos = basePos[i + 2];
			setEndPos();
		}
	}

	return {
		baseTraces,
		...chromData,
	};
}