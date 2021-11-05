#include <mitsuba/core/properties.h>
#include <mitsuba/core/qmc.h>
#include <mitsuba/core/random.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/render/sampler.h>

NAMESPACE_BEGIN(mitsuba)

MTS_VARIANT Sampler<Float, Spectrum>::Sampler(const Properties &props) {
    m_sample_count = props.size_("sample_count", 4);
    m_base_seed    = props.size_("seed", 0);

    m_dimension_index       = 0u;
    m_sample_index          = 0;
    m_samples_per_wavefront = 1;
    m_wavefront_size        = 0;
}

MTS_VARIANT Sampler<Float, Spectrum>::~Sampler() {}

MTS_VARIANT void Sampler<Float, Spectrum>::seed(uint64_t seed_offset) {
    // NotImplementedError("seed");
    // m_wavefront_size  = wavefront_size;
    m_dimension_index = 0u;
    m_sample_index    = 0;
}

MTS_VARIANT void Sampler<Float, Spectrum>::advance() {
    Assert(m_sample_index < (m_sample_count / m_samples_per_wavefront));
    m_dimension_index = 0u;
    m_sample_index++;
}

MTS_VARIANT Float Sampler<Float, Spectrum>::next_1d(Mask) {
    NotImplementedError("next_1d");
}

MTS_VARIANT typename Sampler<Float, Spectrum>::Point2f
Sampler<Float, Spectrum>::next_2d(Mask) {
    NotImplementedError("next_2d");
}

MTS_VARIANT void Sampler<Float, Spectrum>::set_samples_per_wavefront(
    uint32_t samples_per_wavefront) {
    if constexpr (is_scalar_v<Float>)
        Throw("set_samples_per_wavefront should not be used in scalar variants "
              "of the renderer.");

    m_samples_per_wavefront = samples_per_wavefront;
    if (m_sample_count % m_samples_per_wavefront != 0)
        Throw("sample_count should be a multiple of samples_per_wavefront!");
}

MTS_VARIANT typename Sampler<Float, Spectrum>::UInt32
Sampler<Float, Spectrum>::compute_per_sequence_seed(
    uint32_t seed_offset) const {
    UInt32 indices = arange<UInt32>(m_wavefront_size);
    UInt32 sequence_idx =
        m_samples_per_wavefront * (indices / m_samples_per_wavefront);
    return sample_tea_32(UInt32(m_base_seed),
                         sequence_idx + UInt32(seed_offset));
}

MTS_VARIANT typename Sampler<Float, Spectrum>::UInt32
Sampler<Float, Spectrum>::current_sample_index() const {
    // Build an array of offsets for the sample indices in the wavefront
    UInt32 wavefront_sample_offsets = 0;
    if (m_samples_per_wavefront > 1)
        wavefront_sample_offsets =
            arange<UInt32>(m_wavefront_size) % m_samples_per_wavefront;

    return m_sample_index * m_samples_per_wavefront + wavefront_sample_offsets;
}

MTS_VARIANT LowDiscrepancySampler<Float, Spectrum>::LowDiscrepancySampler(
    const Properties &props)
    : Base(props) {
    // Make sure sample_count is power of two and square (e.g. 4, 16, 64,
    // 256, 1024, ...)
    ScalarUInt32 res = 2;
    while (sqr(res) < m_sample_count)
        res = math::round_to_power_of_two(++res);

    // if (m_sample_count != sqr(res))
    //     Log(Warn,
    //         "Sample count should be square and power of two, rounding to
    //         "
    //         "%i",
    //         sqr(res));

    m_sample_count = sqr(res);
}

MTS_VARIANT LowDiscrepancySampler<Float, Spectrum>::LowDiscrepancySampler(
    int res_s, const Properties &props)
    : Base(props) {
    // Make sure sample_count is power of two and square (e.g. 4, 16, 64,
    // 256, 1024, ...)
    m_sample_count = res_s;
    // res            = sqrt(m_sample_count);
    // ScalarUInt32 res = res_s;
    // while (sqr(res) < m_sample_count)
    //     res = math::round_to_power_of_two(++res);

    // if (m_sample_count != sqr(res))
    //     Log(Warn,
    //         "Sample count should be square and power of two, rounding
    // to
    //         "
    //         "%i",
    //         sqr(res));

    // m_sample_count = sqr(res);
}

// MTS_VARIANT ref<Sampler<Float, Spectrum>>
// LowDiscrepancySampler<Float, Spectrum>::clone() override {
//     LowDiscrepancySampler *sampler   = new LowDiscrepancySampler();
//     sampler->m_sample_count          = m_sample_count;
//     sampler->m_samples_per_wavefront = m_samples_per_wavefront;
//     sampler->m_base_seed             = m_base_seed;
//     return sampler;
// }

// void seed(uint64_t seed_offset, size_t wavefront_size) override {
//     Base::seed(seed_offset, wavefront_size);
//     m_scramble_seed = compute_per_sequence_seed(seed_offset);
// }

MTS_VARIANT void
LowDiscrepancySampler<Float, Spectrum>::seed(uint64_t seed_offset) {
    Base::seed(seed_offset);
    m_scramble_seed = compute_per_sequence_seed(seed_offset);
}

MTS_VARIANT Float LowDiscrepancySampler<Float, Spectrum>::next_1d_test() {
    // Assert(seeded());

    uint32_t sample_indices = current_sample_index();
    uint32_t perm_seed      = m_scramble_seed + m_dimension_index++;

    // Shuffle the samples order
    uint32_t i = permute(sample_indices, m_sample_count, perm_seed);

    // Compute scramble value (unique per sequence)
    uint32_t scramble = sample_tea_32(m_scramble_seed, uint32_t(0x48bc48eb));

    return radical_inverse_2(i, scramble);
}

MTS_VARIANT typename LowDiscrepancySampler<Float, Spectrum>::Point2f
LowDiscrepancySampler<Float, Spectrum>::next_2d_test() {
    // Assert(seeded());

    uint32_t sample_indices = current_sample_index();
    uint32_t perm_seed      = m_scramble_seed + m_dimension_index++;

    // Shuffle the samples order
    uint32_t i = permute(sample_indices, m_sample_count, perm_seed);

    // Compute scramble values (unique per sequence) for both axis
    uint32_t scramble_x = sample_tea_32(m_scramble_seed, uint32_t(0x98bc51ab));
    uint32_t scramble_y = sample_tea_32(m_scramble_seed, uint32_t(0x04223e2d));

    Float x = radical_inverse_2(i, scramble_x), y = sobol_2(i, scramble_y);

    return Point2f(x, y);
}

MTS_IMPLEMENT_CLASS_VARIANT(Sampler, Object, "sampler")
MTS_INSTANTIATE_CLASS(Sampler)
MTS_IMPLEMENT_CLASS_VARIANT(LowDiscrepancySampler, Object, "ldsampler")
MTS_INSTANTIATE_CLASS(LowDiscrepancySampler)
NAMESPACE_END(mitsuba)
