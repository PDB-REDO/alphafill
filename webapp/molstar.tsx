import "core-js/stable";
import 'molstar/build/viewer/molstar.css';
import { StructureSelection } from "molstar/lib/mol-model/structure";
import { createPlugin } from 'molstar/lib/mol-plugin-ui';
import { PluginUIContext } from 'molstar/lib/mol-plugin-ui/context';
import { DefaultPluginUISpec } from 'molstar/lib/mol-plugin-ui/spec';
import { Script } from "molstar/lib/mol-script/script";
import { Asset } from 'molstar/lib/mol-util/assets';
import "regenerator-runtime/runtime";

export { PLUGIN_VERSION as version } from 'molstar/lib/mol-plugin/version';
export { setDebugMode, setProductionMode } from 'molstar/lib/mol-util/debug';

export class Viewer {
	plugin: PluginUIContext

	constructor(element: HTMLElement) {
		this.plugin = createPlugin(element, {
			...DefaultPluginUISpec(),
			layout: {
				initial: {
					isExpanded: false,
					showControls: false
				}
			},
			components: {
				remoteState: 'none'
			}
		});
	}

	loadStructureFromString(model: string) {
		return this.plugin.dataTransaction(async () => {
			await this.plugin.clear();

			const data = await this.plugin.builders.data.rawData({ data: model }, { state: { isGhost: true } });
			const trajectory = await this.plugin.builders.structure.parseTrajectory(data, 'mmcif');
	
			await this.plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default', {
				structure: {
					name: 'model',
					params: { }
				},
				showUnitcell: false,
				representationPreset: 'auto'
			});
        });
	}

	loadStructureFromUrl(url: string) {
		return this.plugin.dataTransaction(async () => {
			await this.plugin.clear();

			const data = await this.plugin.builders.data.download({ url: Asset.Url(url), isBinary: false }, { state: { isGhost: true } });
			const trajectory = await this.plugin.builders.structure.parseTrajectory(data, 'mmcif');
	
			await this.plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default', {
				structure: {
					name: 'model',
					params: { }
				},
				showUnitcell: false,
				representationPreset: 'auto'
			});
        });
	}

	handleResize() {
		this.plugin.layout.events.updated.next(void 0);
	}

	selectAsym(asymID: string) {
		const data = this.plugin.managers.structure.hierarchy.current.structures[0].cell.obj.data;
		if (!data) return;

		const selection = Script.getStructureSelection(Q => Q.struct.generator.atomGroups({
			'chain-test': Q.core.rel.eq([asymID, Q.ammp('label_asym_id')])
		}), data);

		if (StructureSelection.isEmpty(selection)) {
			alert(`asym ${asymID} not found!`);
		}
		else {
			const loci = StructureSelection.toLociWithSourceUnits(selection);

			this.plugin.managers.interactivity.lociSelects.select({ loci });
			this.plugin.managers.structure.focus.setFromLoci(loci);
			this.plugin.managers.camera.focusLoci(loci);
		}
	}
}
