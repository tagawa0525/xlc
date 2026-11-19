use std::sync::{Arc, Mutex};

use tracing::field::{Field, Visit};
use tracing::{Event, Level, Subscriber};
use tracing_subscriber::Registry;
use tracing_subscriber::layer::{Context, Layer, SubscriberExt};

/// 受けたイベントの (名前、target、水準、フィールドの「名前=値」の並び)
type Seen = Arc<Mutex<Vec<(String, String, Level, Vec<String>)>>>;

struct Record(Seen);

struct Fields(Vec<String>);

impl Visit for Fields {
    fn record_debug(&mut self, field: &Field, value: &dyn std::fmt::Debug) {
        self.0.push(format!("{}={value:?}", field.name()));
    }
}

impl<S: Subscriber> Layer<S> for Record {
    fn on_event(&self, event: &Event<'_>, _: Context<'_, S>) {
        let meta = event.metadata();
        let mut fields = Fields(Vec::new());
        event.record(&mut fields);
        self.0.lock().unwrap().push((
            meta.name().to_string(),
            meta.target().to_string(),
            *meta.level(),
            fields.0,
        ));
    }
}

fn record(emit: impl FnOnce()) -> Vec<(String, String, Level, Vec<String>)> {
    let seen = Seen::default();
    let subscriber = Registry::default().with(Record(seen.clone()));
    tracing::subscriber::with_default(subscriber, emit);
    let seen = seen.lock().unwrap();
    seen.clone()
}

const NAME: &str = concat!("event ", file!());
const TARGET: &str = module_path!();

#[test]
fn every_level_names_the_event_by_the_calling_file_without_the_line() {
    let seen = record(|| {
        crate::trace!("t");
        crate::debug!("d");
        crate::info!("i");
        crate::warn!("w");
        crate::error!("e");
    });
    let levels: Vec<_> = seen.iter().map(|(_, _, level, _)| *level).collect();
    assert_eq!(
        levels,
        [
            Level::TRACE,
            Level::DEBUG,
            Level::INFO,
            Level::WARN,
            Level::ERROR
        ]
    );
    for (name, target, _, _) in &seen {
        assert_eq!(name, NAME);
        assert_eq!(target, TARGET);
    }
}

#[test]
fn fields_and_messages_take_the_forms_of_the_tracing_macros() {
    let x = 1;
    let text = "a";
    let seen = record(|| {
        crate::info!("message {}", 2);
        crate::info!(x, "with {}", 3);
        crate::info!(k = 4, "keyed");
        crate::info!(?text);
        crate::info!(%text, "shown");
        crate::info!({ x }, "braced {}", 5);
        crate::info!(target: "other", x);
    });
    let fields: Vec<_> = seen
        .iter()
        .map(|(_, _, _, fields)| fields.clone())
        .collect();
    assert_eq!(
        fields,
        [
            vec!["message=message 2".to_string()],
            vec!["message=with 3".to_string(), "x=1".to_string()],
            vec!["message=keyed".to_string(), "k=4".to_string()],
            vec![r#"text="a""#.to_string()],
            vec!["message=shown".to_string(), "text=a".to_string()],
            vec!["message=braced 5".to_string(), "x=1".to_string()],
            vec!["x=1".to_string()],
        ]
    );
    let targets: Vec<_> = seen
        .iter()
        .map(|(_, target, _, _)| target.as_str())
        .collect();
    assert_eq!(
        targets,
        [TARGET, TARGET, TARGET, TARGET, TARGET, TARGET, "other"]
    );
    assert!(seen.iter().all(|(name, _, _, _)| name == NAME));
}
