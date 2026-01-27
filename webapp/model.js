import "core-js/stable";
import 'bootstrap';
import 'pdbe-molstar/build/pdbe-molstar-component';

function updateModel(viewer, cbs, showAllCB) {
	const selected = cbs
		.filter(c => c.checked)
		.map(c => c.getAttribute("data-asym-id"));

	const allChecked = selected.length == cbs.length;

	if (allChecked) {
		showAllCB.indeterminate = false;
		showAllCB.checked = true;
	}
	else if (selected.length == 0) {
		showAllCB.indeterminate = false;
		showAllCB.checked = false;
	}
	else
		showAllCB.indeterminate = true;

	const linkAs = [...document.querySelectorAll('a.optimize-link')];
	linkAs.forEach(a => {
		a.classList.toggle('invisible', selected.length != 1 || a.getAttribute('data-asym-id') != selected[0]);
	});

	return viewer.visual.update({
		customData: {
			url: `${window.location.origin}/v1/aff/${AF_ID}/stripped/${selected.join(',')}/${IDENTITY}`,
			format: "cif",
			binary: false
		},
		bgColor: "white"
	}, true);
}

window.addEventListener('load', () => {

	const showAllCB = document.getElementById('show-all');
	const cbs = [...document.querySelectorAll("tr.transplanted-row input[type='checkbox']")];
	const selected = cbs
		.filter(c => c.checked)
		.map(c => c.getAttribute("data-asym-id"));

	const molstarContainer = document.getElementById("app");
	const viewer = new PDBeMolstarPlugin();

	const options = {
		bgColor: "white",
		customData: {
			url: `${window.location.origin}/v1/aff/${AF_ID}/stripped/${selected.join(',')}/${IDENTITY}`,
			format: "cif",
			binary: false
		},
		hideControls: true
	};

	viewer.render(molstarContainer, options);

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

	const rows = document.querySelectorAll("tr.transplanted-row");
	[...rows].forEach(row => {
		row.addEventListener('click', () => {

			const asymID = row.getAttribute('data-asym-id');
			const cb = row.querySelector("input[type='checkbox']");

			if (cb.checked) {
				viewer.visual.select({ data: [{ struct_asym_id: asymID, color: "#2378de" }] })
					.then(() => viewer.visual.focus([{ struct_asym_id: asymID }]));
			}
			else {
				cb.checked = true;
				updateModel(viewer, cbs, showAllCB)
					.then(() => viewer.visual.focus([{ struct_asym_id: asymID }]))
					.then(() => viewer.visual.select({ data: [{ struct_asym_id: asymID, color: "#2378de" }] }));
			}
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
})
