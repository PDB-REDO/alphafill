import "core-js/stable";
import 'molstar/build/viewer/molstar.css';
import { G3dProvider } from 'molstar/lib/extensions/g3d/format';
import { StructureSelection } from "molstar/lib/mol-model/structure";
import { DownloadStructure } from 'molstar/lib/mol-plugin-state/actions/structure';
import { StructureRepresentationPresetProvider } from 'molstar/lib/mol-plugin-state/builder/structure/representation-preset';
import { DataFormatProvider } from 'molstar/lib/mol-plugin-state/formats/provider';
import { BuiltInTrajectoryFormat } from 'molstar/lib/mol-plugin-state/formats/trajectory';
import { createPlugin } from 'molstar/lib/mol-plugin-ui';
import { PluginUIContext } from 'molstar/lib/mol-plugin-ui/context';
import { DefaultPluginUISpec, PluginUISpec } from 'molstar/lib/mol-plugin-ui/spec';
import { PluginConfig } from 'molstar/lib/mol-plugin/config';
import { PluginLayoutControlsDisplay } from 'molstar/lib/mol-plugin/layout';
import { Script } from "molstar/lib/mol-script/script";
import { Asset } from 'molstar/lib/mol-util/assets';
import "regenerator-runtime/runtime";

export { PLUGIN_VERSION as version } from 'molstar/lib/mol-plugin/version';
export { setDebugMode, setProductionMode } from 'molstar/lib/mol-util/debug';

const DefaultViewerOptions = {
	layoutIsExpanded: true,
	layoutShowControls: true,
	layoutShowRemoteState: true,
	layoutControlsDisplay: 'reactive' as PluginLayoutControlsDisplay,
	layoutShowSequence: true,
	layoutShowLog: true,
	layoutShowLeftPanel: true,
	collapseLeftPanel: false,
	collapseRightPanel: false,
	disableAntialiasing: PluginConfig.General.DisableAntialiasing.defaultValue,
	pixelScale: PluginConfig.General.PixelScale.defaultValue,
	pickScale: PluginConfig.General.PickScale.defaultValue,
	pickPadding: PluginConfig.General.PickPadding.defaultValue,
	enableWboit: PluginConfig.General.EnableWboit.defaultValue,

	viewportShowExpand: PluginConfig.Viewport.ShowExpand.defaultValue,
	viewportShowControls: PluginConfig.Viewport.ShowControls.defaultValue,
	viewportShowSettings: PluginConfig.Viewport.ShowSettings.defaultValue,
	viewportShowSelectionMode: PluginConfig.Viewport.ShowSelectionMode.defaultValue,
	viewportShowAnimation: PluginConfig.Viewport.ShowAnimation.defaultValue,
	pluginStateServer: PluginConfig.State.DefaultServer.defaultValue,
	volumeStreamingServer: PluginConfig.VolumeStreaming.DefaultServer.defaultValue,
	volumeStreamingDisabled: !PluginConfig.VolumeStreaming.Enabled.defaultValue,
	pdbProvider: PluginConfig.Download.DefaultPdbProvider.defaultValue,
	emdbProvider: PluginConfig.Download.DefaultEmdbProvider.defaultValue,
};
type ViewerOptions = typeof DefaultViewerOptions;

export class Viewer {
	plugin: PluginUIContext

	constructor(elementOrId: string | HTMLElement, options: Partial<ViewerOptions> = {}) {
		const o = { ...DefaultViewerOptions, ...options };
		const defaultSpec = DefaultPluginUISpec();

		const spec: PluginUISpec = {
			actions: defaultSpec.actions,
			behaviors: [
				...defaultSpec.behaviors,
			],
			animations: [...defaultSpec.animations || []],
			customParamEditors: defaultSpec.customParamEditors,
			layout: {
				initial: {
					isExpanded: o.layoutIsExpanded,
					showControls: o.layoutShowControls,
					controlsDisplay: o.layoutControlsDisplay,
					regionState: {
						bottom: 'full',
						left: o.collapseLeftPanel ? 'collapsed' : 'full',
						right: o.collapseRightPanel ? 'hidden' : 'full',
						top: 'full',
					}
				},
			},
			components: {
				...defaultSpec.components,
				controls: {
					...defaultSpec.components?.controls,
					top: o.layoutShowSequence ? undefined : 'none',
					bottom: o.layoutShowLog ? undefined : 'none',
					left: o.layoutShowLeftPanel ? undefined : 'none',
				},
				remoteState: o.layoutShowRemoteState ? 'default' : 'none',
			},
			config: [
				[PluginConfig.General.DisableAntialiasing, o.disableAntialiasing],
				[PluginConfig.General.PixelScale, o.pixelScale],
				[PluginConfig.General.PickScale, o.pickScale],
				[PluginConfig.General.PickPadding, o.pickPadding],
				[PluginConfig.General.EnableWboit, o.enableWboit],
				[PluginConfig.Viewport.ShowExpand, o.viewportShowExpand],
				[PluginConfig.Viewport.ShowControls, o.viewportShowControls],
				[PluginConfig.Viewport.ShowSettings, o.viewportShowSettings],
				[PluginConfig.Viewport.ShowSelectionMode, o.viewportShowSelectionMode],
				[PluginConfig.Viewport.ShowAnimation, o.viewportShowAnimation],
				[PluginConfig.State.DefaultServer, o.pluginStateServer],
				[PluginConfig.State.CurrentServer, o.pluginStateServer],
				[PluginConfig.VolumeStreaming.Enabled, false],
			]
		};

		const element = typeof elementOrId === 'string'
			? document.getElementById(elementOrId)
			: elementOrId;
		if (!element) throw new Error(`Could not get element with id '${elementOrId}'`);
		this.plugin = createPlugin(element, spec);
	}

	loadStructureFromUrl(url: string, format: BuiltInTrajectoryFormat = 'mmcif', isBinary = false, options?: LoadStructureOptions) {
		const params = DownloadStructure.createDefaultParams(this.plugin.state.data.root.obj!, this.plugin);
		return this.plugin.runTask(this.plugin.state.data.applyAction(DownloadStructure, {
			source: {
				name: 'url',
				params: {
					url: Asset.Url(url),
					format: format as any,
					isBinary,
					options: { ...params.source.params.options, representationParams: options?.representationParams as any },
				}
			}
		}));
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

export interface LoadStructureOptions {
	representationParams?: StructureRepresentationPresetProvider.CommonParams
}

