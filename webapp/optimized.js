import "core-js/stable";
import "regenerator-runtime/runtime";
import 'bootstrap';
import 'pdbe-molstar/build/pdbe-molstar-component';

window.addEventListener('load', async () => {
	const loadAlert = document.getElementById('load-alert');
	const errorAlert = document.getElementById('error-alert');

	try {
		const r = await fetch(`v1/aff/${AF_ID}/optimized-with-stats/A,${ASYM_ID}`);
		loadAlert.style.display = 'none';
		const data = await r.json();

		if (typeof data.error == 'string')
			throw data.error;

		const molstarContainer = document.getElementById("app");
		const viewer = new PDBeMolstarPlugin();

		viewer.events.loadComplete.subscribe(() => {
			viewer.visual.focus([{ struct_asym_id: ASYM_ID }]);
		});

		const options = {
			bgColor: "white",
			customData: {
				url: `data:${encodeURI(data.model)}`,
				format: "cif",
				binary: false
			},
			hideControls: true
		};

		viewer.render(molstarContainer, options);

		const formatter = new Intl.NumberFormat('en-US', {
			minimumFractionDigits: 2,
			maximumFractionDigits: 2,
		});

		document.getElementById("clash-original").textContent = formatter.format(+result.clash.before);
		document.getElementById("clash-optimized").textContent = formatter.format(+result.clash.after);

		const link = document.getElementById('model-link');
		link.href = "data:text/plain;charset=utf-8," + encodeURIComponent(data.model);
		link.download = `${AF_ID}-${ASYM_ID}-optimized.cif`;

		document.getElementById('link-table').classList.remove('invisible');

		document.getElementById('model').classList.remove('invisible');
	} catch (err) {
		loadAlert.style.display = 'none';
		errorAlert.style.display = '';
		const msg = document.getElementById('error-message');
		msg.textContent = err;
	}
})
