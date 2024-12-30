import numpy as np
import pynwb
from spikeinterface import BaseRecording, BaseSorting
from spikeextractors import RecordingExtractor, SortingExtractor
from hdmf.data_utils import DataChunkIterator
from hdmf.backends.hdf5.h5_utils import H5DataIO



def add_electrical_series(
    recording: SpikeInterfaceRecording,
    nwbfile=None,
    metadata: dict = None,
    segment_index: int = 0,
    starting_time: Optional[float] = None,
    use_times: bool = False,
    write_as: str = "raw",
    es_key: str = None,
    write_scaled: bool = False,
    compression: Optional[str] = "gzip",
    compression_opts: Optional[int] = None,
    iterator_type: Optional[str] = None,
    iterator_opts: Optional[dict] = None,
):
    """
    Auxiliary static method for nwbextractor.
    Adds traces from recording object as ElectricalSeries to nwbfile object.
    Parameters
    ----------
    recording: SpikeInterfaceRecording
    nwbfile: NWBFile
        nwb file to which the recording information is to be added
    metadata: dict
        metadata info for constructing the nwb file (optional).
        Should be of the format
            metadata['Ecephys']['ElectricalSeries'] = dict(
                name=my_name,
                description=my_description
            )
    segment_index : int
        The recording segment to add to the NWBFile.
    starting_time: float (optional)
        Sets the starting time of the ElectricalSeries to a manually set value.
        Increments timestamps if use_times is True.
    use_times: bool (optional, defaults to False)
        If True, the times are saved to the nwb file using recording.frame_to_time(). If False (defualut),
        the sampling rate is used.
    write_as: str (optional, defaults to 'raw')
        How to save the traces data in the nwb file. Options:
        - 'raw' will save it in acquisition
        - 'processed' will save it as FilteredEphys, in a processing module
        - 'lfp' will save it as LFP, in a processing module
    es_key: str (optional)
        Key in metadata dictionary containing metadata info for the specific electrical series
    write_scaled: bool (optional, defaults to True)
        If True, writes the scaled traces (return_scaled=True)
    compression: str (optional, defaults to "gzip")
        Type of compression to use. Valid types are "gzip" and "lzf".
        Set to None to disable all compression.
    compression_opts: int (optional, defaults to 4)
        Only applies to compression="gzip". Controls the level of the GZIP.
    iterator_type: str (optional, defaults to 'v2')
        The type of DataChunkIterator to use.
        'v1' is the original DataChunkIterator of the hdmf data_utils.
        'v2' is the locally developed RecordingExtractorDataChunkIterator, which offers full control over chunking.
    iterator_opts: dict (optional)
        Dictionary of options for the RecordingExtractorDataChunkIterator (iterator_type='v2').
        Valid options are
            buffer_gb: float (optional, defaults to 1 GB)
                Recommended to be as much free RAM as available. Automatically calculates suitable buffer shape.
            chunk_mb: float (optional, defaults to 1 MB)
                Should be below 1 MB. Automatically calculates suitable chunk shape.
        If manual specification of buffer_shape and chunk_shape are desired, these may be specified as well.
    Missing keys in an element of metadata['Ecephys']['ElectrodeGroup'] will be auto-populated with defaults
    whenever possible.
    """
    if isinstance(recording, RecordingExtractor):
        checked_recording = OldToNewRecording(oldapi_recording_extractor=recording)
    else:
        checked_recording = recording
    if nwbfile is not None:
        assert isinstance(nwbfile, pynwb.NWBFile), "'nwbfile' should be of type pynwb.NWBFile!"
    assert compression is None or compression in [
        "gzip",
        "lzf",
    ], "Invalid compression type ({compression})! Choose one of 'gzip', 'lzf', or None."

    if not nwbfile.electrodes:
        add_electrodes(recording, nwbfile, metadata)
    assert write_as in [
        "raw",
        "processed",
        "lfp",
    ], f"'write_as' should be 'raw', 'processed' or 'lfp', but instead received value {write_as}"

    if compression == "gzip":
        if compression_opts is None:
            compression_opts = 4
        else:
            assert compression_opts in range(
                10
            ), "compression type is 'gzip', but specified compression_opts is not an integer between 0 and 9!"
    elif compression == "lzf" and compression_opts is not None:
        warn(f"compression_opts ({compression_opts}) were passed, but compression type is 'lzf'! Ignoring options.")
        compression_opts = None
    if iterator_opts is None:
        iterator_opts = dict()
    if write_as == "raw":
        eseries_kwargs = dict(
            name="ElectricalSeries_raw",
            description="Raw acquired data",
            comments="Generated from SpikeInterface::NwbRecordingExtractor",
        )
    elif write_as == "processed":
        eseries_kwargs = dict(
            name="ElectricalSeries_processed",
            description="Processed data",
            comments="Generated from SpikeInterface::NwbRecordingExtractor",
        )
        ecephys_mod = get_module(
            nwbfile=nwbfile,
            name="ecephys",
            description="Intermediate data from extracellular electrophysiology recordings, e.g., LFP.",
        )
        if "Processed" not in ecephys_mod.data_interfaces:
            ecephys_mod.add(pynwb.ecephys.FilteredEphys(name="Processed"))
    elif write_as == "lfp":
        eseries_kwargs = dict(
            name="ElectricalSeries_lfp",
            description="Processed data - LFP",
            comments="Generated from SpikeInterface::NwbRecordingExtractor",
        )
        ecephys_mod = get_module(
            nwbfile=nwbfile,
            name="ecephys",
            description="Intermediate data from extracellular electrophysiology recordings, e.g., LFP.",
        )
        if "LFP" not in ecephys_mod.data_interfaces:
            ecephys_mod.add(pynwb.ecephys.LFP(name="LFP"))
    if metadata is not None and "Ecephys" in metadata and es_key is not None:
        assert es_key in metadata["Ecephys"], f"metadata['Ecephys'] dictionary does not contain key '{es_key}'"
        eseries_kwargs.update(metadata["Ecephys"][es_key])
    if write_as == "raw":
        assert (
            eseries_kwargs["name"] not in nwbfile.acquisition
        ), f"Raw ElectricalSeries '{eseries_kwargs['name']}' is already written in the NWBFile!"
    elif write_as == "processed":
        assert (
            eseries_kwargs["name"] not in nwbfile.processing["ecephys"].data_interfaces["Processed"].electrical_series
        ), f"Processed ElectricalSeries '{eseries_kwargs['name']}' is already written in the NWBFile!"
    elif write_as == "lfp":
        assert (
            eseries_kwargs["name"] not in nwbfile.processing["ecephys"].data_interfaces["LFP"].electrical_series
        ), f"LFP ElectricalSeries '{eseries_kwargs['name']}' is already written in the NWBFile!"

    # Indexes by channel ids if they are integer or by indices otherwise.
    channel_name_array = checked_recording.get_channel_ids()
    if np.issubdtype(channel_name_array.dtype, np.integer):
        channel_indices = channel_name_array
    else:
        channel_indices = checked_recording.ids_to_indices(channel_name_array)

    table_ids = [list(nwbfile.electrodes.id[:]).index(id) for id in channel_indices]

    electrode_table_region = nwbfile.create_electrode_table_region(
        region=table_ids, description="electrode_table_region"
    )
    eseries_kwargs.update(electrodes=electrode_table_region)

    # channels gains - for RecordingExtractor, these are values to cast traces to uV.
    # For nwb, the conversions (gains) cast the data to Volts.
    # To get traces in Volts we take data*channel_conversion*conversion.
    channel_conversion = checked_recording.get_channel_gains()
    channel_offset = checked_recording.get_channel_offsets()
    if write_scaled or channel_conversion is None:
        eseries_kwargs.update(conversion=1e-6)
    else:
        if len(np.unique(channel_conversion)) == 1:  # if all gains are equal
            eseries_kwargs.update(conversion=channel_conversion[0] * 1e-6)
        else:
            eseries_kwargs.update(conversion=1e-6)
            eseries_kwargs.update(channel_conversion=channel_conversion)
    if iterator_type is None or iterator_type == "v2":
        ephys_data = SpikeInterfaceRecordingDataChunkIterator(
            recording=checked_recording,
            segment_index=segment_index,
            return_scaled=write_scaled,
            **iterator_opts,
        )
    elif iterator_type == "v1":
        if isinstance(checked_recording.get_traces(end_frame=5, return_scaled=write_scaled), np.memmap) and np.all(
            channel_offset == 0
        ):
            ephys_data = DataChunkIterator(
                data=checked_recording.get_traces(return_scaled=write_scaled), **iterator_opts
            )
        else:
            raise ValueError("iterator_type='v1' only supports memmapable trace types! Use iterator_type='v2' instead.")
    else:
        raise NotImplementedError(f"iterator_type ({iterator_type}) should be either 'v1' or 'v2' (recommended)!")
    eseries_kwargs.update(data=H5DataIO(data=ephys_data, compression=compression, compression_opts=compression_opts))

    if not use_times and starting_time is None:
        eseries_kwargs.update(starting_time=float(checked_recording.get_times(segment_index=segment_index)[0]))
    elif not use_times and starting_time is not None:
        eseries_kwargs.update(starting_time=starting_time)
    if not use_times:
        eseries_kwargs.update(rate=float(recording.get_sampling_frequency()))
    elif not use_times and starting_time is not None:
        eseries_kwargs.update(rate=float(checked_recording.get_sampling_frequency()))
    elif use_times and starting_time is not None:
        eseries_kwargs.update(
            timestamps=H5DataIO(
                data=starting_time
                + checked_recording.get_times()[
                    np.arange(checked_recording.get_num_samples(segment_index=segment_index))
                ],
                compression=compression,
                compression_opts=compression_opts,
            )
        )
    elif use_times and starting_time is None:
        eseries_kwargs.update(
            timestamps=H5DataIO(
                data=checked_recording.get_times()[
                    np.arange(checked_recording.get_num_samples(segment_index=segment_index))
                ],
                compression=compression,
                compression_opts=compression_opts,
            )
        )
    es = pynwb.ecephys.ElectricalSeries(**eseries_kwargs)
    if write_as == "raw":
        nwbfile.add_acquisition(es)
    elif write_as == "processed":
        ecephys_mod.data_interfaces["Processed"].add_electrical_series(es)
    elif write_as == "lfp":
        ecephys_mod.data_interfaces["LFP"].add_electrical_series(es)
        
