import "core-js/stable";
import "regenerator-runtime/runtime";

import 'bootstrap';

import { Viewer } from "./molstar";

function updateModel(viewer, cbs, showAllCB) {
	const selected = cbs
		.filter(c => c.checked)
		.map(c => c.getAttribute("data-asym-id"));
	
	const allChecked = selected.length == cbs.length;

	if (allChecked)
	{
		showAllCB.indeterminate = false;
		showAllCB.checked = true;
	}
	else if (selected.length == 0)
	{
		showAllCB.indeterminate = false;
		showAllCB.checked = false;
	}
	else
		showAllCB.indeterminate = true;
	
	selected.push('A');
	
	return viewer.loadStructureFromUrl(`v1/aff/${AF_ID}/stripped/${selected.join(',')}/${IDENTITY}`);
}

window.addEventListener('load', () => {

	const viewer = new Viewer(document.getElementById('app'));

	const showAllCB = document.getElementById('show-all');

	const cbs = [...document.querySelectorAll("tr.transplanted-row input[type='checkbox']")];

	showAllCB.addEventListener('change', () => {
		const checked = showAllCB.checked;

		cbs.forEach(cb => cb.checked = checked);

		updateModel(viewer, cbs, showAllCB);
	});

	cbs.forEach(cb => {
		cb.addEventListener('click', (evt) => evt.stopPropagation());

		cb.addEventListener('change', () => updateModel(viewer, cbs, showAllCB));
	});

	// viewer.loadStructureFromUrl(`/v1/aff/${AF_ID}`)
	updateModel(viewer, cbs, showAllCB)
		.then(() => {
			const rows = document.querySelectorAll("tr.transplanted-row");
			[...rows].forEach(row => {
				row.addEventListener('click', () => {

					const asymID = row.getAttribute('data-asym-id');
					const cb = row.querySelector("input[type='checkbox']");

					if (cb.checked)
						viewer.selectAsym(asymID);
					else {
						cb.checked = true;
						updateModel(viewer, cbs, showAllCB).then(() => viewer.selectAsym(asymID));
					}
				});
			});
		});

	// identity buttons

	const ibs = [...document.querySelectorAll("input[type='radio']")];
	ibs.forEach(ib => {
		const identity = ib.getAttribute('data-identity');
		ib.addEventListener('click', () => window.location = `model?id=${AF_ID}&identity=${identity}`)
	});

})
