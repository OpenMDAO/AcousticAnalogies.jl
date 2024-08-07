const CompactSourceElement = CompactF1ASourceElement
Base.@deprecate CompactSourceElement CompactF1ASourceElement

Base.@deprecate f1a(se::CompactF1ASourceElement, obs::AbstractAcousticObserver, t_obs) noise(se::CompactF1ASourceElement, obs::AbstractAcousticObserver, t_obs)
Base.@deprecate f1a(se::CompactF1ASourceElement, obs::AbstractAcousticObserver) noise(se::CompactF1ASourceElement, obs::AbstractAcousticObserver)

Base.@deprecate source_elements_ccblade(rotor, sections, ops, outputs, area_per_chord2, period, num_src_times, positive_x_rotation) f1a_source_elements_ccblade(rotor, sections, ops, outputs, area_per_chord2, period, num_src_times, positive_x_rotation)
