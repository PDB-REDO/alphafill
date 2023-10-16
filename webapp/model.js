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

	const linkAs = [...document.querySelectorAll('a.optimize-link')];
	linkAs.forEach(a => {
		a.classList.toggle('invisible', selected.length != 1 || a.getAttribute('data-asym-id') != selected[0]);
	});

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


	const links = [...document.querySelectorAll("tr.transplanted-row a")];
	links.forEach(link => {
		link.addEventListener('click', (evt) => evt.stopPropagation());
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


	// download button
	const downloadBtn = document.getElementById('structure-with-selected-ligands');
	downloadBtn.addEventListener('click', (e) => {
		e.preventDefault();

		const selected = cbs
			.filter(c => c.checked)
			.map(c => c.getAttribute("data-asym-id"));
	
		window.location = `v1/aff/${AF_ID}/stripped/${selected.join(',')}/${IDENTITY}`;
	});

	// // Update PAE matrices
	// const cv = [...document.querySelectorAll("td.pae-matrix canvas")];
	// cv.forEach(c => {
	// 	const ctx = c.getContext("2d");

	// 	const pae = JSON.parse(c.getAttribute("data-pae"));
		
	// 	const dim = Math.sqrt(pae.length);
	// 	const scale = 32.0 / dim;

	// 	ctx.scale(scale, scale);

	// 	let i = 0;
	// 	for (let x = 0; x < dim; ++x)
	// 	{
	// 		for (let y = 0; y < dim; ++y)
	// 		{
	// 			const v = 1.0 * pae[i];
	// 			i += 1;

	// 			console.log(v);
				
	// 			const clr = `hsl(220,60%,${Math.round(25 + 50 * (v / 32.0))}%)`;

	// 			ctx.fillStyle = clr;
	// 			ctx.fillRect(x, y, 1, 1);
	// 		}
	// 	}
	// });


})
